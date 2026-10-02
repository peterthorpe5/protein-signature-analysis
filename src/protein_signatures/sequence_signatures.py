"""Amino-acid feature generation for protein-signature discovery."""

from __future__ import annotations

import logging
from collections import Counter
from heapq import nsmallest

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
    vocabulary_policy: str = "strict",
    maximum_candidates: int | None = None,
) -> tuple[FeatureRecord, ...]:
    """Generate protein-level presence features for amino-acid k-mers.

    Args:
        sequences: Authoritative protein records.
        discovery_protein_ids: Proteins eligible to define the feature vocabulary.
            Retained discovery-defined k-mers are subsequently scanned across every
            supplied sequence, including held-out validation proteins.
        lengths: Positive k-mer lengths.
        minimum_proteins: Minimum discovery-protein prevalence retained.
        maximum_features: Raw candidate limit in strict mode; maximum number of
            retained features in prevalence-ranked mode.
        vocabulary_policy: ``strict`` retains all eligible candidates or fails;
            ``prevalence_ranked`` shares the feature budget across lengths,
            then fills unused slots with the most prevalent discovery strings.
            Equal-prevalence strings are ordered lexically.
        maximum_candidates: Raw candidate safeguard required for ranked mode.

    Returns:
        Deterministically ordered protein/k-mer presence records.

    Raises:
        InputValidationError: If settings are invalid or the safeguard is exceeded.
    """
    if not lengths or any(length < 1 for length in lengths):
        raise InputValidationError("At least one positive k-mer length is required.")
    if len(lengths) != len(set(lengths)):
        raise InputValidationError("K-mer lengths must not contain duplicates.")
    if minimum_proteins < 1 or maximum_features < 1:
        raise InputValidationError("K-mer prevalence and feature limits must be positive.")
    if not isinstance(vocabulary_policy, str) or vocabulary_policy not in {
        "strict",
        "prevalence_ranked",
    }:
        raise InputValidationError("Unknown k-mer vocabulary policy.")
    if vocabulary_policy == "prevalence_ranked":
        if (
            isinstance(maximum_candidates, bool)
            or not isinstance(maximum_candidates, int)
            or maximum_candidates < maximum_features
        ):
            raise InputValidationError(
                "Ranked k-mer selection requires maximum_candidates >= maximum_features."
            )
    elif maximum_candidates is not None:
        raise InputValidationError("maximum_candidates applies only to ranked k-mer selection.")
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

    prevalence: Counter[str] = Counter()
    candidate_limit = (
        maximum_features if vocabulary_policy == "strict" else maximum_candidates
    )
    for record in sequences:
        if record.protein_id not in discovery_protein_ids:
            continue
        observed = frozenset(
            f"k{length}:{kmer}"
            for length in lengths
            for kmer in sequence_kmers(sequence=record.sequence, length=length)
        )
        prevalence.update(observed)
        if len(prevalence) > candidate_limit:
            setting = (
                "analysis.maximum_kmer_features"
                if vocabulary_policy == "strict"
                else "analysis.maximum_kmer_candidates"
            )
            raise InputValidationError(
                "The unique discovery k-mer safeguard was exceeded "
                f"({len(prevalence):,} > {candidate_limit:,}); reduce k-mer lengths or "
                f"increase {setting} deliberately."
            )
    if vocabulary_policy == "strict":
        retained = frozenset(
            kmer for kmer, count in prevalence.items() if count >= minimum_proteins
        )
        eligible_count = len(retained)
    else:
        selected: set[str] = set()
        eligible_count = 0
        ordered_lengths = sorted(lengths)
        equal_budget, extra_slots = divmod(maximum_features, len(ordered_lengths))
        counts_by_length: dict[int, int] = {}
        for index, length in enumerate(ordered_lengths):
            prefix = f"k{length}:"
            count = sum(
                observed >= minimum_proteins and kmer.startswith(prefix)
                for kmer, observed in prevalence.items()
            )
            counts_by_length[length] = count
            eligible_count += count
            budget = equal_budget + (index >= len(ordered_lengths) - extra_slots)
            if budget:
                selected.update(
                    nsmallest(
                        budget,
                        (
                            kmer
                            for kmer, observed in prevalence.items()
                            if observed >= minimum_proteins and kmer.startswith(prefix)
                        ),
                        key=lambda kmer: (-prevalence[kmer], kmer),
                    )
                )
        remaining = maximum_features - len(selected)
        if remaining and eligible_count > len(selected):
            selected.update(
                nsmallest(
                    remaining,
                    (
                        kmer
                        for kmer, count in prevalence.items()
                        if count >= minimum_proteins and kmer not in selected
                    ),
                    key=lambda kmer: (-prevalence[kmer], kmer),
                )
            )
        for length in ordered_lengths:
            prefix = f"k{length}:"
            LOGGER.info(
                "K-mer vocabulary length=%d eligible=%d retained=%d",
                length,
                counts_by_length[length],
                sum(kmer.startswith(prefix) for kmer in selected),
            )
        retained = frozenset(selected)
        if eligible_count > len(retained):
            LOGGER.warning(
                "Prevalence-ranked discovery vocabulary retained %d/%d eligible k-mers "
                "from %d raw candidates; unselected strings cannot be tested",
                len(retained),
                eligible_count,
                len(prevalence),
            )
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
