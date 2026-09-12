"""Leakage-resistant provenance contracts for generic protein features."""

from __future__ import annotations

import logging
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from types import MappingProxyType

from .checksums import sha256_text
from .errors import InputValidationError
from .models import FeatureRecord, SequenceRecord

LOGGER = logging.getLogger(__name__)

ASSESSED_WITH_FEATURE = "ASSESSED_WITH_FEATURE"
ASSESSED_NO_FEATURE = "ASSESSED_NO_FEATURE"
FEATURE_ASSESSMENT_STATUSES = frozenset(
    {
        ASSESSED_WITH_FEATURE,
        ASSESSED_NO_FEATURE,
        "NOT_ASSESSED",
        "FAILED",
        "EXCLUDED",
    }
)
FIXED_EXTERNAL = "FIXED_EXTERNAL"
DISCOVERY_DERIVED = "DISCOVERY_DERIVED"
DISCOVERY_SUPERVISED = "DISCOVERY_SUPERVISED"
ALL_DATA_EXPLORATORY = "ALL_DATA_EXPLORATORY"
FEATURE_DERIVATION_SCOPES = frozenset(
    {
        FIXED_EXTERNAL,
        DISCOVERY_DERIVED,
        DISCOVERY_SUPERVISED,
        ALL_DATA_EXPLORATORY,
    }
)


@dataclass(frozen=True)
class ImportedFeatureEvidence:
    """Prepared generic evidence safe for confirmatory inference.

    Attributes:
        confirmatory_features: Positive rows whose definitions did not use held-out
            campaign data.
        assessment_records: Every imported row, including assessed absence and
            unknown/failure states, retained for the published audit ledger.
        assessment_universes: Successfully assessed proteins for each confirmatory
            feature key.
        exploratory_feature_keys: Supervised discovery or all-data-derived feature
            definitions excluded from confirmatory association and modelling.
    """

    confirmatory_features: tuple[FeatureRecord, ...]
    assessment_records: tuple[FeatureRecord, ...]
    assessment_universes: Mapping[tuple[str, str], frozenset[str]]
    exploratory_feature_keys: frozenset[tuple[str, str]]


def sequence_cohort_sha256(
    *, sequences: tuple[SequenceRecord, ...], protein_ids: Collection[str]
) -> str:
    """Hash an exact, order-independent protein-sequence cohort.

    Args:
        sequences: Authoritative protein inventory.
        protein_ids: Non-empty subset included in the derivation cohort.

    Returns:
        SHA-256 of canonical sorted ``protein_id``, sequence-digest rows.

    Raises:
        InputValidationError: If records are duplicated or cohort IDs are invalid.
    """
    if isinstance(protein_ids, (str, bytes)) or not isinstance(protein_ids, Collection):
        raise InputValidationError("protein_ids must be a collection of protein identifiers.")
    selected = frozenset(protein_ids)
    if not selected or any(not isinstance(item, str) or not item.strip() for item in selected):
        raise InputValidationError("A cohort requires one or more non-empty protein IDs.")
    by_id = {record.protein_id: record for record in sequences}
    if len(by_id) != len(sequences):
        raise InputValidationError("Protein identifiers must be unique for cohort hashing.")
    unknown = sorted(selected - frozenset(by_id))
    if unknown:
        raise InputValidationError(
            f"Cohort protein identifiers are absent from the supplied sequences: {unknown[:10]}"
        )
    canonical = "".join(
        f"{protein_id}\t{by_id[protein_id].sequence_sha256}\n" for protein_id in sorted(selected)
    )
    return sha256_text(text=canonical)


def prepare_imported_feature_evidence(
    *,
    records: tuple[FeatureRecord, ...],
    sequences: tuple[SequenceRecord, ...],
    discovery_protein_ids: Collection[str],
) -> ImportedFeatureEvidence:
    """Validate derivation cohorts and isolate confirmatory generic features.

    A ``FIXED_EXTERNAL`` definition is prespecified independently of this campaign.
    A ``DISCOVERY_DERIVED`` definition attests that derivation was label-blind and
    must name the exact discovery-cohort digest. A ``DISCOVERY_SUPERVISED``
    definition may use discovery target/background labels and must name that same
    digest; it is retained for audit but excluded from v0.1 inference.
    ``ALL_DATA_EXPLORATORY`` definitions must name the full campaign cohort and are
    likewise audit-only. Successfully assessed presence and absence rows define
    feature-specific denominators for confirmatory definitions; missing, failed,
    excluded and not-assessed rows remain unknown. The package verifies declared
    cohort identity but cannot inspect whether an external producer used labels.

    Args:
        records: Strictly parsed generic feature-assessment rows.
        sequences: Authoritative protein inventory.
        discovery_protein_ids: Frozen discovery set used by feature learning.

    Returns:
        Positive confirmatory features, the complete ledger and assessed universes.

    Raises:
        InputValidationError: If derivation metadata or cohort digests disagree.
    """
    discovery_digest = sequence_cohort_sha256(
        sequences=sequences, protein_ids=discovery_protein_ids
    )
    full_digest = sequence_cohort_sha256(
        sequences=sequences,
        protein_ids=frozenset(record.protein_id for record in sequences),
    )
    metadata_by_key: dict[tuple[str, str], tuple[str, str, str]] = {}
    universes: dict[tuple[str, str], set[str]] = {}
    positives: list[FeatureRecord] = []
    exploratory: set[tuple[str, str]] = set()
    for record in records:
        key = (record.feature_type, record.feature_id)
        metadata = (
            record.derivation_scope,
            record.feature_definition_sha256,
            record.derivation_cohort_sha256,
        )
        previous = metadata_by_key.setdefault(key, metadata)
        if previous != metadata:
            raise InputValidationError(
                f"Conflicting derivation metadata for generic feature {key!r}."
            )
        if record.derivation_scope == FIXED_EXTERNAL:
            if record.derivation_cohort_sha256:
                raise InputValidationError(
                    f"FIXED_EXTERNAL feature {key!r} must leave derivation_cohort_sha256 blank."
                )
        elif record.derivation_scope == DISCOVERY_DERIVED:
            if record.derivation_cohort_sha256 != discovery_digest:
                raise InputValidationError(
                    f"DISCOVERY_DERIVED feature {key!r} does not match the exact "
                    "discovery cohort digest."
                )
        elif record.derivation_scope == DISCOVERY_SUPERVISED:
            if record.derivation_cohort_sha256 != discovery_digest:
                raise InputValidationError(
                    f"DISCOVERY_SUPERVISED feature {key!r} does not match the exact "
                    "discovery cohort digest."
                )
            exploratory.add(key)
            continue
        elif record.derivation_scope == ALL_DATA_EXPLORATORY:
            if record.derivation_cohort_sha256 != full_digest:
                raise InputValidationError(
                    f"ALL_DATA_EXPLORATORY feature {key!r} does not match the exact "
                    "full campaign cohort digest."
                )
            exploratory.add(key)
            continue
        else:
            raise InputValidationError(
                f"Unknown derivation scope for generic feature {key!r}: "
                f"{record.derivation_scope!r}."
            )
        universes.setdefault(key, set())
        if record.evidence_status in {ASSESSED_WITH_FEATURE, ASSESSED_NO_FEATURE}:
            universes[key].add(record.protein_id)
        if record.evidence_status == ASSESSED_WITH_FEATURE:
            positives.append(record)
    frozen_universes = {
        key: frozenset(protein_ids) for key, protein_ids in sorted(universes.items())
    }
    LOGGER.info(
        "Prepared %d confirmatory generic feature rows across %d feature definitions; "
        "excluded %d supervised-discovery or all-data definitions from inference",
        len(positives),
        len(frozen_universes),
        len(exploratory),
    )
    return ImportedFeatureEvidence(
        confirmatory_features=tuple(positives),
        assessment_records=records,
        assessment_universes=MappingProxyType(frozen_universes),
        exploratory_feature_keys=frozenset(exploratory),
    )
