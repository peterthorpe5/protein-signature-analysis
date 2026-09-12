"""Auditable AlphaFold Protein Structure Database model acquisition."""

from __future__ import annotations

import json
import logging
import os
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from .checksums import sha256_file
from .errors import InputValidationError
from .io_utils import iter_tsv
from .models import (
    AlphaFoldAcquisition,
    AlphaFoldRequest,
    AlphaFoldSettings,
    FoldEvidenceStatus,
    SequenceRecord,
    StructureAnalysisEligibility,
    StructureRecord,
)
from .validation import validate_identifier

LOGGER = logging.getLogger(__name__)
ALPHAFOLD_API_ROOT = "https://alphafold.ebi.ac.uk/api/prediction"
ACCESSION_FIELDS = ("protein_id", "uniprot_accession")
_USER_AGENT = "protein-signature-analysis/0.1.0"
_MAX_METADATA_BYTES = 5 * 1024 * 1024
_MAX_MODEL_BYTES = 512 * 1024 * 1024


def read_alphafold_requests(
    *, path: Path, protein_ids: frozenset[str]
) -> tuple[AlphaFoldRequest, ...]:
    """Read explicit campaign-protein to UniProt accession mappings.

    Args:
        path: Two-column AlphaFold accession TSV.
        protein_ids: Authoritative FASTA identifiers.

    Returns:
        Unique, deterministic acquisition requests.

    Raises:
        InputValidationError: If a protein is unknown or mapped more than once.
    """

    requests: list[AlphaFoldRequest] = []
    observed: set[str] = set()
    for row in iter_tsv(path=path, required_fields=ACCESSION_FIELDS):
        protein_id = validate_identifier(value=row["protein_id"], field_name="protein_id")
        if protein_id not in protein_ids:
            raise InputValidationError(
                f"AlphaFold accession mapping references unknown protein: {protein_id!r}"
            )
        if protein_id in observed:
            raise InputValidationError(
                f"Protein has multiple AlphaFold accession mappings: {protein_id!r}"
            )
        observed.add(protein_id)
        accession = validate_identifier(
            value=row["uniprot_accession"], field_name="uniprot_accession"
        ).upper()
        requests.append(AlphaFoldRequest(protein_id=protein_id, uniprot_accession=accession))
    return tuple(sorted(requests, key=lambda item: item.protein_id))


def acquire_alphafold_models(
    *,
    requests: tuple[AlphaFoldRequest, ...],
    sequences: tuple[SequenceRecord, ...],
    settings: AlphaFoldSettings,
) -> tuple[tuple[AlphaFoldAcquisition, ...], tuple[StructureRecord, ...]]:
    """Acquire sequence-verified AlphaFold DB PDB models into a checksum cache.

    Args:
        requests: Explicit UniProt accession mappings.
        sequences: Authoritative campaign sequences.
        settings: Cache, timeout, retry and quality settings.

    Returns:
        Acquisition outcomes and successfully acquired structure records.

    Raises:
        InputValidationError: If acquisition is enabled without requests.
    """

    if not settings.enabled:
        outcomes = tuple(
            _outcome(
                request=request,
                status="NOT_SELECTED",
                api_url=_metadata_url(accession=request.uniprot_accession),
                message="AlphaFold acquisition is disabled in the campaign configuration.",
            )
            for request in requests
        )
        return outcomes, ()
    if not requests:
        raise InputValidationError(
            "AlphaFold acquisition is enabled but inputs.alphafold_accessions is absent or empty."
        )
    settings.cache_dir.mkdir(parents=True, exist_ok=True)
    sequences_by_id = {item.protein_id: item for item in sequences}
    outcomes: list[AlphaFoldAcquisition] = []
    structures: list[StructureRecord] = []
    for request in requests:
        outcome, structure = _acquire_one(
            request=request,
            sequence=sequences_by_id[request.protein_id],
            settings=settings,
        )
        outcomes.append(outcome)
        if structure is not None:
            structures.append(structure)
    return (
        tuple(sorted(outcomes, key=lambda item: item.protein_id)),
        tuple(sorted(structures, key=lambda item: item.structure_id)),
    )


def parse_mean_plddt_from_pdb(*, path: Path) -> float | None:
    """Calculate mean AlphaFold pLDDT from PDB C-alpha B-factor fields.

    Args:
        path: AlphaFold-format PDB coordinate file.

    Returns:
        Mean C-alpha pLDDT, or ``None`` when no valid atoms are present.
    """

    values: list[float] = []
    try:
        with Path(path).open(mode="r", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("ENDMDL"):
                    break
                if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                    continue
                alternate = line[16:17]
                if alternate not in {" ", "A"}:
                    continue
                try:
                    value = float(line[60:66].strip())
                except ValueError:
                    continue
                if 0.0 <= value <= 100.0:
                    values.append(value)
    except (OSError, UnicodeError) as error:
        raise InputValidationError(f"Could not parse AlphaFold PDB {path}: {error}") from error
    return sum(values) / len(values) if values else None


def _acquire_one(
    *, request: AlphaFoldRequest, sequence: SequenceRecord, settings: AlphaFoldSettings
) -> tuple[AlphaFoldAcquisition, StructureRecord | None]:
    """Acquire one AlphaFold Database prediction without aborting the campaign.

    Args:
        request: Protein/accession mapping.
        sequence: Authoritative campaign sequence.
        settings: Acquisition settings.

    Returns:
        Acquisition outcome and optional retained structure. The structure can be
        ineligible for analysis while remaining available for provenance review.
    """

    api_url = _metadata_url(accession=request.uniprot_accession)
    try:
        metadata_bytes = _request_bytes(
            url=api_url,
            timeout_seconds=settings.timeout_seconds,
            retries=settings.retries,
            maximum_bytes=_MAX_METADATA_BYTES,
            allowed_host_suffix="ebi.ac.uk",
        )
        metadata = _select_metadata(payload=metadata_bytes, accession=request.uniprot_accession)
        structure_id = validate_identifier(
            value=metadata.get("entryId"), field_name="AlphaFold entryId"
        )
        version = str(metadata.get("latestVersion", "")).strip()
        remote_sequence = str(metadata.get("uniprotSequence", "")).replace("\n", "").upper()
        sequence_match = remote_sequence == sequence.sequence if remote_sequence else None
        if sequence_match is False:
            return (
                _outcome(
                    request=request,
                    status="SEQUENCE_MISMATCH",
                    api_url=api_url,
                    structure_id=structure_id,
                    model_version=version,
                    sequence_match=False,
                    message="AlphaFold DB sequence differs from the authoritative campaign FASTA.",
                ),
                None,
            )
        pdb_url = _validated_download_url(value=metadata.get("pdbUrl"))
        destination = settings.cache_dir / f"{structure_id}-model_v{version or 'unknown'}.pdb"
        if not destination.is_file() or destination.stat().st_size == 0:
            payload = _request_bytes(
                url=pdb_url,
                timeout_seconds=settings.timeout_seconds,
                retries=settings.retries,
                maximum_bytes=_MAX_MODEL_BYTES,
                allowed_host_suffix="ebi.ac.uk",
            )
            _write_bytes_atomic(path=destination, payload=payload)
        digest = sha256_file(path=destination)
        mean_plddt = parse_mean_plddt_from_pdb(path=destination)
        if sequence_match is None:
            status = "ACQUIRED_SEQUENCE_UNVERIFIED"
            eligibility = StructureAnalysisEligibility.INELIGIBLE_SEQUENCE_UNVERIFIED
        elif mean_plddt is None:
            status = "ACQUIRED_CONFIDENCE_UNAVAILABLE"
            eligibility = StructureAnalysisEligibility.INELIGIBLE_CONFIDENCE_UNAVAILABLE
        elif mean_plddt < settings.minimum_mean_plddt:
            status = "ACQUIRED_LOW_CONFIDENCE"
            eligibility = StructureAnalysisEligibility.INELIGIBLE_LOW_CONFIDENCE
        else:
            status = "ACQUIRED"
            eligibility = StructureAnalysisEligibility.ELIGIBLE
        outcome = _outcome(
            request=request,
            status=status,
            api_url=api_url,
            structure_id=structure_id,
            model_version=version,
            coordinate_path=destination,
            coordinate_sha256=digest,
            sequence_match=sequence_match,
            mean_plddt=mean_plddt,
            message="AlphaFold DB model cached with content checksum.",
        )
        structure = StructureRecord(
            protein_id=request.protein_id,
            structure_id=structure_id,
            structure_source="AlphaFoldDB",
            structure_version=version,
            coordinate_path=destination,
            coordinate_sha256=digest,
            availability_status="AVAILABLE",
            mean_confidence=mean_plddt,
            fold_id="",
            fold_name="",
            fold_authority="",
            fold_authority_version="",
            fold_evidence_reference="",
            fold_evidence_status=FoldEvidenceStatus.NOT_ASSESSED,
            analysis_eligibility_status=eligibility,
            comparison_universe_ids=(),
        )
        LOGGER.info(
            "Retained AlphaFold DB model %s for %s with eligibility %s",
            structure_id,
            request.protein_id,
            eligibility.value,
        )
        return outcome, structure
    except urllib.error.HTTPError as error:
        status = "MODEL_NOT_AVAILABLE" if error.code == 404 else "FAILED"
        message = f"AlphaFold DB HTTP response {error.code}."
    except (
        InputValidationError,
        OSError,
        UnicodeError,
        ValueError,
        urllib.error.URLError,
    ) as error:
        status = "FAILED"
        message = f"AlphaFold acquisition failed: {error}"
    LOGGER.warning("%s Protein %s (%s)", message, request.protein_id, request.uniprot_accession)
    return _outcome(request=request, status=status, api_url=api_url, message=message), None


def _request_bytes(
    *,
    url: str,
    timeout_seconds: float,
    retries: int,
    maximum_bytes: int = _MAX_MODEL_BYTES,
    allowed_host_suffix: str | None = None,
) -> bytes:
    """Retrieve one HTTPS resource with bounded retries.

    Args:
        url: Validated HTTPS URL.
        timeout_seconds: Per-request timeout.
        retries: Number of retries after the first attempt.
        maximum_bytes: Hard response-body ceiling.
        allowed_host_suffix: Optional required hostname or parent-domain suffix,
            checked again after redirects.

    Returns:
        Non-empty response bytes.

    Raises:
        OSError: If all attempts fail or the response is empty.
        urllib.error.HTTPError: For a final HTTP failure.
    """

    if maximum_bytes < 1:
        raise InputValidationError("AlphaFold response byte ceiling must be positive.")
    if allowed_host_suffix is not None:
        _require_allowed_host(url=url, allowed_host_suffix=allowed_host_suffix)
    request = urllib.request.Request(url=url, headers={"User-Agent": _USER_AGENT})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                final_url = str(response.geturl())
                if allowed_host_suffix is not None:
                    _require_allowed_host(
                        url=final_url,
                        allowed_host_suffix=allowed_host_suffix,
                    )
                content_length = response.headers.get("Content-Length")
                if content_length is not None:
                    try:
                        declared_size = int(content_length)
                    except ValueError as error:
                        raise OSError(
                            f"Invalid Content-Length from {final_url}: {content_length!r}"
                        ) from error
                    if declared_size > maximum_bytes:
                        raise OSError(
                            f"Response from {final_url} exceeds the {maximum_bytes}-byte ceiling."
                        )
                payload = response.read(maximum_bytes + 1)
            if not payload:
                raise OSError(f"Empty response from {url}")
            if len(payload) > maximum_bytes:
                raise OSError(f"Response from {url} exceeds the {maximum_bytes}-byte ceiling.")
            return payload
        except urllib.error.HTTPError as error:
            if error.code < 500 or attempt == retries:
                raise
        except (OSError, urllib.error.URLError):
            if attempt == retries:
                raise
        time.sleep(min(2**attempt, 8))
    raise OSError(f"Request failed without a response: {url}")


def _require_allowed_host(*, url: str, allowed_host_suffix: str) -> None:
    """Require an HTTPS URL within one exact host/domain suffix.

    Args:
        url: URL before or after redirection.
        allowed_host_suffix: Allowed host or parent-domain suffix.

    Raises:
        InputValidationError: If the URL is not HTTPS or leaves the allowed domain.
    """

    parsed = urllib.parse.urlparse(url)
    hostname = (parsed.hostname or "").casefold()
    allowed = allowed_host_suffix.casefold().lstrip(".")
    if parsed.scheme != "https" or not (hostname == allowed or hostname.endswith(f".{allowed}")):
        raise InputValidationError(f"Unexpected AlphaFold DB response URL: {url!r}")


def _select_metadata(*, payload: bytes, accession: str) -> dict[str, Any]:
    """Select the latest matching model from AlphaFold DB API metadata.

    Args:
        payload: JSON response bytes.
        accession: Requested UniProt accession.

    Returns:
        Selected metadata mapping.

    Raises:
        InputValidationError: If metadata is malformed or lacks a match.
    """

    try:
        document = json.loads(payload.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise InputValidationError(f"Invalid AlphaFold DB metadata JSON: {error}") from error
    if not isinstance(document, list):
        raise InputValidationError("AlphaFold DB metadata response must be a list.")
    candidates = [
        item
        for item in document
        if isinstance(item, dict)
        and str(item.get("uniprotAccession", "")).upper() == accession.upper()
    ]
    if not candidates:
        raise InputValidationError(f"AlphaFold DB metadata has no model for {accession!r}.")
    return max(candidates, key=lambda item: int(item.get("latestVersion", 0)))


def _validated_download_url(*, value: Any) -> str:
    """Validate an AlphaFold DB model download URL.

    Args:
        value: Candidate URL from official API metadata.

    Returns:
        Validated HTTPS URL.

    Raises:
        InputValidationError: If scheme or host is unexpected.
    """

    url = str(value or "").strip()
    try:
        _require_allowed_host(url=url, allowed_host_suffix="ebi.ac.uk")
    except InputValidationError as error:
        raise InputValidationError(f"Unexpected AlphaFold DB download URL: {url!r}") from error
    return url


def _metadata_url(*, accession: str) -> str:
    """Build the official AlphaFold DB metadata endpoint for an accession.

    Args:
        accession: Validated UniProt accession.

    Returns:
        HTTPS API URL.
    """

    safe_accession = urllib.parse.quote(accession, safe="")
    return f"{ALPHAFOLD_API_ROOT}/{safe_accession}"


def _write_bytes_atomic(*, path: Path, payload: bytes) -> None:
    """Write non-empty model bytes atomically into the cache.

    Args:
        path: Final cache path.
        payload: Downloaded coordinate bytes.

    Raises:
        InputValidationError: If the payload is empty.
    """

    if not payload:
        raise InputValidationError(f"Refusing to cache an empty coordinate model at {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, mode="wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except OSError:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _outcome(
    *,
    request: AlphaFoldRequest,
    status: str,
    api_url: str,
    message: str,
    structure_id: str = "",
    model_version: str = "",
    coordinate_path: Path | None = None,
    coordinate_sha256: str = "",
    sequence_match: bool | None = None,
    mean_plddt: float | None = None,
) -> AlphaFoldAcquisition:
    """Construct a consistently shaped AlphaFold acquisition outcome.

    Args:
        request: Protein/accession mapping.
        status: Controlled acquisition status.
        api_url: Metadata endpoint.
        message: Human-readable diagnostic.
        structure_id: Acquired model identifier.
        model_version: AlphaFold DB model version.
        coordinate_path: Cached PDB path.
        coordinate_sha256: Cached PDB checksum.
        sequence_match: FASTA/API sequence comparison result.
        mean_plddt: Mean model confidence.

    Returns:
        Immutable acquisition record.
    """

    return AlphaFoldAcquisition(
        protein_id=request.protein_id,
        uniprot_accession=request.uniprot_accession,
        acquisition_status=status,
        structure_id=structure_id,
        model_version=model_version,
        api_url=api_url,
        coordinate_path=coordinate_path,
        coordinate_sha256=coordinate_sha256,
        sequence_match=sequence_match,
        mean_plddt=mean_plddt,
        message=message,
    )
