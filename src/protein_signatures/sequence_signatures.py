"""Amino-acid feature generation for protein-signature discovery."""

from __future__ import annotations

import logging
from collections import Counter

from .checksums import sha256_text
from .errors import InputValidationError
from .fasta import sequence_kmers
from .feature_provenance import DISCOVERY_DERIVED, sequence_cohort_sha256
from .models import FeatureRecord, SequenceRecord

LOGGER = logging.getLogger(__name__)


def build_kmer_features(
    *,
    sequences: tuple[SequenceRecord, ...],
    discovery_protein_ids: frozenset[str],
    lengths: tuple[int, ...],
    minimum_proteins: int,
    maximum_features: int,
) -> tuple[FeatureRecord, ...]:
    """Generate protein-level presence features for amino-acid k-mers.

    Args:
        sequences: Authoritative protein records.
        discovery_protein_ids: Proteins eligible to define the feature vocabulary.
            Retained discovery-defined k-mers are subsequently scanned across every
            supplied sequence, including held-out validation proteins.
        lengths: Positive k-mer lengths.
        minimum_proteins: Minimum discovery-protein prevalence retained.
        maximum_features: Maximum unique discovery candidate k-mers before failing
            safely.

    Returns:
        Deterministically ordered protein/k-mer presence records.

    Raises:
        InputValidationError: If settings are invalid or the safeguard is exceeded.
    """
    if not lengths or any(length < 1 for length in lengths):
        raise InputValidationError("At least one positive k-mer length is required.")
    if minimum_proteins < 1 or maximum_features < 1:
        raise InputValidationError("K-mer prevalence and feature limits must be positive.")
    sequence_ids = [record.protein_id for record in sequences]
    if len(set(sequence_ids)) != len(sequence_ids):
        raise InputValidationError("Protein identifiers must be unique for k-mer generation.")
    if not discovery_protein_ids:
        raise InputValidationError(
            "At least one discovery protein is required to define the k-mer vocabulary."
        )
    unknown_discovery_ids = sorted(discovery_protein_ids - frozenset(sequence_ids))
    if unknown_discovery_ids:
        preview = ", ".join(repr(item) for item in unknown_discovery_ids[:5])
        suffix = " ..." if len(unknown_discovery_ids) > 5 else ""
        raise InputValidationError(
            "Discovery protein identifiers are absent from the supplied sequences: "
            f"{preview}{suffix}"
        )
    derivation_cohort_sha256 = sequence_cohort_sha256(
        sequences=sequences,
        protein_ids=discovery_protein_ids,
    )

    discovery_kmers: dict[str, frozenset[str]] = {}
    prevalence: Counter[str] = Counter()
    for record in sequences:
        if record.protein_id not in discovery_protein_ids:
            continue
        observed = frozenset(
            f"k{length}:{kmer}"
            for length in lengths
            for kmer in sequence_kmers(sequence=record.sequence, length=length)
        )
        discovery_kmers[record.protein_id] = observed
        prevalence.update(observed)
        if len(prevalence) > maximum_features:
            raise InputValidationError(
                "The unique discovery k-mer safeguard was exceeded "
                f"({len(prevalence):,} > {maximum_features:,}); reduce k-mer lengths or "
                "increase analysis.maximum_kmer_features deliberately."
            )
    retained = frozenset(kmer for kmer, count in prevalence.items() if count >= minimum_proteins)
    definition_digests = {
        kmer: sha256_text(
            text=(
                "protein-signature-analysis feature definition v1\n"
                "feature_type=AMINO_ACID_KMER\n"
                f"feature_id={kmer}\n"
                "algorithm=exact_protein_level_presence\n"
            )
        )
        for kmer in retained
    }
    features: list[FeatureRecord] = []
    for record in sorted(sequences, key=lambda item: item.protein_id):
        observed = discovery_kmers.get(record.protein_id)
        if observed is None:
            observed = frozenset(
                f"k{length}:{kmer}"
                for length in lengths
                for kmer in sequence_kmers(sequence=record.sequence, length=length)
            )
        for kmer in sorted(observed & retained):
            features.append(
                FeatureRecord(
                    protein_id=record.protein_id,
                    feature_type="AMINO_ACID_KMER",
                    feature_id=kmer,
                    feature_name=kmer.replace(":", " ", 1),
                    start=None,
                    end=None,
                    evidence_status="DERIVED",
                    evidence_source="protein-signature-analysis",
                    evidence_reference="discovery_defined_exact_protein_level_kmer_presence",
                    derivation_scope=DISCOVERY_DERIVED,
                    feature_definition_sha256=definition_digests[kmer],
                    derivation_cohort_sha256=derivation_cohort_sha256,
                )
            )
    LOGGER.info(
        "Retained %d discovery-defined amino-acid k-mers represented by %d "
        "protein-feature rows across all partitions",
        len(retained),
        len(features),
    )
    return tuple(features)
