"""Optional ESMFold and annotation services with exact-sequence provenance."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
import time
from pathlib import Path
from urllib.parse import urlencode

import requests

from protein_signatures.errors import ExternalToolError, InputValidationError
from protein_signatures.io_utils import write_json_atomic, write_tsv_atomic

from .analysis import iter_fasta, run_command
from .catalogue import identifiers
from .config import WorkbenchConfig, https_endpoint
from .structures import MODEL_BYTES, model_residues

LOGGER = logging.getLogger(__name__)
INTERPRO = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
ANNOTATION_FIELDS = (
    "sequence_id",
    "sequence_sha256",
    "authority",
    "accession",
    "description",
    "start",
    "end",
    "score",
    "interpro_accession",
    "go_terms",
    "evidence_type",
)
AUDIT_FIELDS = ("sequence_id", "sequence_sha256", "status", "hit_count")


def select_proteins(
    *, source: Path, selected_ids: list[str], maximum: int, scan_limit: int
) -> dict[str, str]:
    """Select a complete exact set from a protein FASTA without sampling.

    Args:
        source: Unaligned protein FASTA.
        selected_ids: Exact selected identifiers; blank selects the whole input.
        maximum: Maximum selected proteins.
        scan_limit: Maximum records to inspect.

    Returns:
        Selected exact proteins in FASTA order.
    """
    wanted = set(identifiers(values=selected_ids)) if selected_ids else set()
    if len(wanted) > maximum:
        raise InputValidationError(f"Select at most {maximum} exact protein IDs.")
    proteins = {}
    for count, (identifier, sequence) in enumerate(iter_fasta(path=source), start=1):
        if count > scan_limit:
            raise InputValidationError("Input FASTA exceeds the workbench scan bound.")
        if wanted and identifier not in wanted:
            continue
        if identifier in proteins:
            raise InputValidationError(f"Duplicate selected protein ID: {identifier}")
        proteins[identifier] = sequence
        if len(proteins) > maximum:
            raise InputValidationError(f"Select at most {maximum} exact protein IDs.")
    if wanted - proteins.keys():
        raise InputValidationError(
            "At least one exact selected protein ID is absent from the FASTA."
        )
    return proteins


def provider_request(
    *,
    url: str,
    method: str = "GET",
    body: bytes | None = None,
    headers: dict[str, str] | None = None,
    maximum: int = MODEL_BYTES,
) -> bytes:
    """Use one trusted HTTPS endpoint with bounded response size and time.

    Args:
        url: Operator-owned or fixed provider URL.
        method: GET or POST.
        body: Explicit request body.
        headers: Content type and optional server-held credentials.
        maximum: Maximum response size in bytes.

    Returns:
        Complete response bytes.

    Raises:
        ExternalToolError: If transport, service status or response bounds fail.
    """
    https_endpoint(value=url)
    deadline = time.monotonic() + 150
    try:
        with requests.request(
            method=method,
            url=url,
            data=body,
            headers=headers,
            timeout=(10, 120),
            allow_redirects=False,
            stream=True,
        ) as response:
            if not 200 <= response.status_code < 300:
                raise ExternalToolError(f"Protein provider returned HTTP {response.status_code}.")
            result = bytearray()
            for chunk in response.iter_content(chunk_size=65536):
                result.extend(chunk)
                if len(result) > maximum or time.monotonic() > deadline:
                    raise ExternalToolError(
                        "Protein provider response exceeded its time or size bound."
                    )
    except requests.RequestException:
        raise ExternalToolError("Could not reach the configured protein provider.") from None
    return bytes(result)


def folding_fasta(*, proteins: dict[str, str], maximum_length: int) -> bytes:
    """Encode the exact bounded ESMFold input set shared by local and portable requests."""
    if not 1 <= len(proteins) <= 20 or any(
        len(sequence) > maximum_length or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX]+", sequence)
        for sequence in proteins.values()
    ):
        raise InputValidationError(
            f"ESMFold accepts 1–20 canonical proteins (or X) up to {maximum_length} residues."
        )
    return "".join(
        f">{identifier}\n{sequence}\n" for identifier, sequence in proteins.items()
    ).encode()


def predict_models(
    *, config: WorkbenchConfig, source: Path, parameters: dict, output: Path
) -> None:
    """Submit a small, explicitly approved set to a configured ESMFold endpoint.

    Args:
        config: Trusted API settings and limits.
        source: Exact protein FASTA input.
        parameters: Selected IDs and explicit sequence submission consent.
        output: Owned job directory.
    """
    if not config.esm_endpoint or parameters.get("remote_allowed") is not True:
        raise InputValidationError(
            "ESMFold needs a configured endpoint and explicit sequence submission."
        )
    if config.esm_format not in {"sequence", "json_inputs"}:
        raise InputValidationError("Unsupported ESMFold API input format.")
    proteins = select_proteins(
        source=source,
        selected_ids=parameters.get("selected_ids", []),
        maximum=20,
        scan_limit=config.max_sequences,
    )
    fasta = folding_fasta(proteins=proteins, maximum_length=config.esm_max_length)
    token = os.environ.get(config.esm_token_env) if config.esm_token_env else None
    if config.esm_token_env and not token:
        raise InputValidationError("The worker is missing the configured ESMFold token variable.")
    content_type = "application/json" if config.esm_format == "json_inputs" else "text/plain"
    headers = {"Content-Type": content_type}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    models = output / "models"
    models.mkdir()
    rows, confidence = [], []
    (output / "model_inputs.fasta").write_bytes(data=fasta)
    for index, (identifier, sequence) in enumerate(proteins.items(), start=1):
        body = (
            json.dumps({"inputs": sequence}).encode()
            if content_type == "application/json"
            else (sequence.encode())
        )
        pdb = provider_request(url=config.esm_endpoint, method="POST", body=body, headers=headers)
        if pdb.lstrip().startswith((b"{", b"[")):
            try:
                decoded = json.loads(pdb)
                if isinstance(decoded, list) and len(decoded) == 1:
                    decoded = decoded[0]
                pdb = decoded["pdb"].encode()
            except (AttributeError, KeyError, TypeError, ValueError) as error:
                raise ExternalToolError("ESMFold response did not include a PDB model.") from error
        suffix, observed, records = model_residues(content=pdb, maximum=config.esm_max_length)
        if (
            suffix != ".pdb"
            or observed != sequence
            or any(not 0 <= row["bfactor"] <= 100 for row in records)
        ):
            raise InputValidationError(
                "Predicted PDB residues or confidence do not match the input."
            )
        name = f"models/model_{index:04d}.pdb"
        (output / name).write_bytes(data=pdb)
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        rows.append(
            {
                "sequence_id": identifier,
                "sequence_sha256": digest,
                "model_file": name,
                "model_sha256": hashlib.sha256(pdb).hexdigest(),
                "model_length": len(observed),
                "mean_reported_confidence": sum(row["bfactor"] for row in records) / len(records),
            }
        )
        confidence.extend({"sequence_id": identifier, **row} for row in records)
    write_tsv_atomic(path=output / "models.tsv", records=rows, fieldnames=tuple(rows[0]))
    write_tsv_atomic(
        path=output / "model_confidence.tsv", records=confidence, fieldnames=tuple(confidence[0])
    )
    write_json_atomic(
        path=output / "esmfold_provenance.json",
        value={
            "endpoint": config.esm_endpoint,
            "input_format": config.esm_format,
            "provider_version": "SERVICE_VERSION_UNVERIFIED",
            "confidence": "Provider-reported PDB B-factor; no experimental validation implied",
        },
    )
    LOGGER.info("Saved %d exact-sequence ESMFold models", len(rows))


def interpro_rows(*, content: str, proteins: dict[str, str]) -> list[dict]:
    """Validate InterPro TSV sequence identities and inclusive positions.

    Args:
        content: Provider result.
        proteins: Exact submitted sequences by identifier.

    Returns:
        Normalised positional predictions.
    """
    rows = []
    for fields in csv.reader(content.splitlines(), delimiter="\t"):
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < 11 or fields[0] not in proteins:
            raise InputValidationError("InterPro returned a malformed row or an unknown protein.")
        identifier, sequence = fields[0], proteins[fields[0]]
        if fields[1].lower() != hashlib.md5(
            sequence.encode(), usedforsecurity=False
        ).hexdigest() or fields[2] != str(len(sequence)):
            raise InputValidationError("InterPro result length or sequence checksum differs.")
        try:
            start, end = int(fields[6]), int(fields[7])
        except ValueError as error:
            raise InputValidationError("InterPro feature positions must be integers.") from error
        if not 1 <= start <= end <= len(sequence):
            raise InputValidationError("InterPro feature lies outside the submitted sequence.")
        rows.append(
            {
                "sequence_id": identifier,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "authority": fields[3],
                "accession": fields[4],
                "description": fields[5],
                "start": start,
                "end": end,
                "score": fields[8],
                "interpro_accession": fields[11] if len(fields) > 11 else "",
                "go_terms": fields[13] if len(fields) > 13 else "",
                "evidence_type": "SEQUENCE_SIGNATURE_PREDICTION",
            }
        )
    return rows


def eggnog_rows(*, path: Path, proteins: dict[str, str]) -> list[dict]:
    """Read eggNOG orthology-derived evidence without inventing feature coordinates.

    Args:
        path: emapper annotations table.
        proteins: Exact submitted sequences.

    Returns:
        Non-positional predicted functional evidence.
    """
    rows, header = [], None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#query\t"):
                header = line[1:].rstrip("\r\n").split("\t")
            elif line.strip() and not line.startswith("#"):
                fields = line.rstrip("\r\n").split("\t")
                if not header or len(header) != len(fields) or fields[0] not in proteins:
                    raise InputValidationError("eggNOG returned an invalid header or protein ID.")
                row = dict(zip(header, fields, strict=True))
                identifier = fields[0]
                rows.append(
                    {
                        "sequence_id": identifier,
                        "sequence_sha256": hashlib.sha256(
                            proteins[identifier].encode()
                        ).hexdigest(),
                        "authority": "eggNOG",
                        "accession": row.get("seed_ortholog", ""),
                        "description": row.get("Description", ""),
                        "score": row.get("score", ""),
                        "go_terms": row.get("GOs", ""),
                        "evidence_type": "ORTHOLOGY_TRANSFER_PREDICTION",
                    }
                )
    if header is None:
        raise InputValidationError("eggNOG output does not contain its annotation header.")
    return rows


def remote_interpro(
    *, config: WorkbenchConfig, proteins: dict[str, str], output: Path
) -> list[dict]:
    """Run bounded EBI jobs with saved service IDs for troubleshooting.

    Args:
        config: Contact email and job time limit.
        proteins: Exact transport-safe protein identifiers and sequences.
        output: Owned job directory for service provenance.

    Returns:
        Validated predicted features from all submitted proteins.
    """
    if not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", config.interpro_email):
        raise InputValidationError("Configure an interpro_email before using the EBI service.")
    rows = []
    deadline = time.monotonic() + min(config.command_timeout, 7200)
    for identifier, sequence in proteins.items():
        if time.monotonic() >= deadline:
            raise ExternalToolError("InterProScan polling time expired; see the saved service IDs.")
        body = urlencode(
            {
                "email": config.interpro_email,
                "sequence": f">{identifier}\n{sequence}\n",
                "goterms": "true",
            }
        ).encode()
        service_id = (
            provider_request(
                url=f"{INTERPRO}/run",
                method="POST",
                body=body,
                maximum=1024,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            .decode()
            .strip()
        )
        if not re.fullmatch(r"[A-Za-z0-9_-]+", service_id):
            raise ExternalToolError("InterProScan returned an invalid service identifier.")
        with (output / "interpro_service_jobs.jsonl").open(mode="a", encoding="utf-8") as audit:
            audit.write(json.dumps({"id": identifier, "service_id": service_id}) + "\n")
        while time.monotonic() < deadline:
            status = (
                provider_request(url=f"{INTERPRO}/status/{service_id}", maximum=1024)
                .decode()
                .strip()
            )
            if status == "FINISHED":
                payload = provider_request(url=f"{INTERPRO}/result/{service_id}/tsv")
                (output / f"interpro_{identifier}.tsv").write_bytes(data=payload)
                rows.extend(
                    interpro_rows(content=payload.decode(), proteins={identifier: sequence})
                )
                break
            if status not in {"PENDING", "RUNNING", "QUEUED"}:
                raise ExternalToolError(f"InterProScan service job {service_id} failed.")
            time.sleep(min(10, max(0, deadline - time.monotonic())))
        else:
            raise ExternalToolError("InterProScan polling time expired; see the saved service IDs.")
    return rows


def annotate_proteins(
    *, config: WorkbenchConfig, source: Path, parameters: dict, output: Path
) -> None:
    """Run InterProScan or eggNOG and retain explicit per-protein assessment states.

    Args:
        config: Provider settings, database paths and resource limits.
        source: Unaligned proteins.
        parameters: Selected IDs, provider, taxonomy and remote consent.
        output: New job output.
    """
    engine = parameters.get("engine", "interpro_local")
    if engine not in {"interpro_local", "interpro_remote", "eggnog"}:
        raise InputValidationError("Choose a supported annotation provider.")
    if engine == "interpro_remote" and parameters.get("remote_allowed") is not True:
        raise InputValidationError("EBI annotation requires explicit sequence submission consent.")
    proteins = select_proteins(
        source=source,
        selected_ids=parameters.get("selected_ids", []),
        maximum=10 if engine == "interpro_remote" else 500,
        scan_limit=config.max_sequences,
    )
    if any(len(sequence) > 10000 for sequence in proteins.values()):
        raise InputValidationError("Annotation inputs must be at most 10,000 residues per protein.")
    aliases = {f"p{index:06d}": name for index, name in enumerate(proteins, start=1)}
    safe_proteins = {key: proteins[name] for key, name in aliases.items()}
    normalised = output / "annotation_input.fasta"
    normalised.write_text(
        data="".join(f">{key}\n{seq}\n" for key, seq in proteins.items()), encoding="utf-8"
    )
    transport = output / "annotation_tool_input.fasta"
    transport.write_text(
        data="".join(f">{key}\n{seq}\n" for key, seq in safe_proteins.items()), encoding="utf-8"
    )
    write_tsv_atomic(
        path=output / "annotation_id_map.tsv",
        fieldnames=("tool_id", "sequence_id"),
        records=[{"tool_id": key, "sequence_id": name} for key, name in aliases.items()],
    )
    if engine == "interpro_remote":
        rows = remote_interpro(config=config, proteins=safe_proteins, output=output)
    elif engine == "interpro_local":
        run_command(
            config=config,
            output=output,
            arguments=[
                config.interproscan,
                "-i",
                str(transport),
                "-b",
                str(output / "interpro"),
                "-f",
                "TSV",
                "-cpu",
                str(config.threads),
                "-goterms",
                "-pa",
            ],
        )
        raw = output / "interpro.tsv"
        if raw.stat().st_size > 100 * 1024 * 1024:
            raise InputValidationError("InterPro result exceeds the 100 MiB parse bound.")
        rows = interpro_rows(content=raw.read_text(encoding="utf-8"), proteins=safe_proteins)
    else:
        scope, taxon = parameters.get("tax_scope", "auto"), parameters.get("target_taxon", "")
        if (
            not isinstance(scope, str)
            or not re.fullmatch(r"[A-Za-z0-9_]+", scope)
            or not isinstance(taxon, str)
            or (taxon and not re.fullmatch(r"[1-9][0-9]{0,9}", taxon))
            or not config.eggnog_data
            or not Path(config.eggnog_data).is_dir()
        ):
            raise InputValidationError("Configure eggnog_data and select a valid taxonomy scope.")
        args = [
            config.emapper,
            "-i",
            str(transport),
            "--itype",
            "proteins",
            "-m",
            "diamond",
            "--output",
            "annotation",
            "--output_dir",
            str(output),
            "--data_dir",
            config.eggnog_data,
            "--cpu",
            str(config.threads),
            "--tax_scope",
            scope,
            "--go_evidence",
            "non-electronic",
        ]
        if taxon:
            args.extend(("--target_taxa", taxon))
        run_command(config=config, output=output, arguments=args)
        rows = eggnog_rows(path=output / "annotation.emapper.annotations", proteins=safe_proteins)
    counts = dict.fromkeys(proteins, 0)
    for row in rows:
        row["sequence_id"] = aliases[row["sequence_id"]]
        counts[row["sequence_id"]] += 1
    write_tsv_atomic(
        path=output / "annotations.tsv",
        fieldnames=ANNOTATION_FIELDS,
        records=({name: row.get(name) for name in ANNOTATION_FIELDS} for row in rows),
    )
    write_tsv_atomic(
        path=output / "annotation_audit.tsv",
        fieldnames=AUDIT_FIELDS,
        records=[
            {
                "sequence_id": name,
                "sequence_sha256": hashlib.sha256(proteins[name].encode()).hexdigest(),
                "status": "ASSESSED_WITH_HIT" if count else "ASSESSED_NO_HIT",
                "hit_count": count,
            }
            for name, count in counts.items()
        ],
    )
    write_json_atomic(
        path=output / "annotation_provenance.json",
        value={
            "engine": engine,
            "parameters": parameters,
            "interpretation": "Predicted evidence; no-hit does not establish absence of function.",
        },
    )
    LOGGER.info(
        "Annotated %d proteins with %s (%d prediction rows)", len(proteins), engine, len(rows)
    )
