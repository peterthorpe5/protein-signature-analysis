"""Validation helpers for feature-specific assessment universes.

Feature rows record observed presence, but absence is only meaningful when the
corresponding assay or annotation process successfully assessed that protein.
This module validates the explicit per-feature universes used to distinguish an
assessed absence from an unknown value.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Collection, Mapping
from typing import TypeAlias

from .errors import InputValidationError
from .feature_provenance import ALL_DATA_EXPLORATORY, DISCOVERY_SUPERVISED
from .models import (
    DomainAssessment,
    DomainAssessmentStatus,
    FeatureRecord,
)

LOGGER = logging.getLogger(__name__)

FeatureKey: TypeAlias = tuple[str, str]
FeatureAssessmentUniverses: TypeAlias = Mapping[FeatureKey, Collection[str]]
NormalisedAssessmentUniverses: TypeAlias = dict[FeatureKey, frozenset[str]]

TECHNICAL_FEATURE_TYPES = frozenset(
    {
        "EVIDENCE_AVAILABILITY",
        "STRUCTURE_AVAILABLE",
    }
)


def normalise_feature_keys(*, keys: Collection[FeatureKey], context: str) -> frozenset[FeatureKey]:
    """Validate an explicit collection of feature keys.

    Args:
        keys: Candidate feature-key collection.
        context: Human-readable name used in validation errors.

    Returns:
        Frozen validated feature keys.

    Raises:
        InputValidationError: If the collection or a key is malformed.
    """

    if isinstance(keys, (str, bytes)) or not isinstance(keys, Collection):
        raise InputValidationError(f"{context.capitalize()} must be a feature-key collection.")
    return frozenset(_validated_feature_key(value=key) for key in keys)


def validate_feature_inference_scopes(
    *,
    features: Collection[FeatureRecord],
    exploratory_feature_keys: Collection[FeatureKey],
    context: str,
) -> frozenset[FeatureKey]:
    """Fail closed on nonconfirmatory feature-derivation declarations.

    Blank derivation scopes on native typed evidence are not interpreted here.
    An explicit all-data scope must also be named in the caller's exploratory
    exclusion/classification set. Discovery-supervised definitions are audit-only
    in v0.1 and cannot enter either association or explainable modelling.

    Args:
        features: Feature records proposed for an inference stage.
        exploratory_feature_keys: Explicit all-data-derived feature keys.
        context: Human-readable stage name used in validation errors.

    Returns:
        Validated exploratory feature keys.

    Raises:
        InputValidationError: If keys are unknown, an all-data declaration is not
            protected explicitly, or a supervised-discovery definition is supplied.
    """

    exploratory = normalise_feature_keys(
        keys=exploratory_feature_keys,
        context=f"{context} exploratory feature keys",
    )
    feature_keys = {(feature.feature_type, feature.feature_id) for feature in features}
    unknown = exploratory - feature_keys
    if unknown:
        raise InputValidationError(
            f"{context.capitalize()} exploratory feature keys do not occur in the "
            f"feature records: {sorted(unknown)[:10]}"
        )
    supervised = {
        (feature.feature_type, feature.feature_id)
        for feature in features
        if feature.derivation_scope == DISCOVERY_SUPERVISED
    }
    if supervised:
        raise InputValidationError(
            f"{context.capitalize()} cannot use DISCOVERY_SUPERVISED audit-only "
            f"feature definitions in v0.1: {sorted(supervised)[:10]}"
        )
    declared_all_data = {
        (feature.feature_type, feature.feature_id)
        for feature in features
        if feature.derivation_scope == ALL_DATA_EXPLORATORY
    }
    unprotected = declared_all_data - exploratory
    if unprotected:
        raise InputValidationError(
            f"{context.capitalize()} requires every ALL_DATA_EXPLORATORY feature to be "
            f"declared exploratory: {sorted(unprotected)[:10]}"
        )
    return exploratory


def compose_feature_assessment_universes(
    *,
    protein_ids: Collection[str],
    features: Collection[FeatureRecord],
    domain_assessments: Collection[DomainAssessment] = (),
    explicit_assessment_universes: Collection[FeatureAssessmentUniverses] = (),
    universally_assessed_feature_types: Collection[str] = ("AMINO_ACID_KMER",),
) -> NormalisedAssessmentUniverses:
    """Compose safe feature-specific universes from campaign evidence ledgers.

    The composer never infers a negative merely because a positive feature row
    is absent. Native sequence feature types explicitly named as universal use
    one shared immutable FASTA inventory rather than per-feature copies. Domain
    features use successful authority assessments; caller-supplied universes
    (for example from imported feature assessments or
    ``derive_structure_assessment_universes``) are merged before strict
    validation. This composer never expands a structural universe from retained
    pairwise hits: an explicit completed-search universe is required upstream.

    Args:
        protein_ids: Authoritative campaign protein identifiers.
        features: Final positive feature records entering analysis.
        domain_assessments: Protein/authority domain assessment ledger.
        explicit_assessment_universes: Zero or more explicit per-feature maps,
            such as imported-feature assessment states and structural universes
            derived from structure records plus complete Foldseek membership.
        universally_assessed_feature_types: Feature types whose derivation
            inspected every authoritative protein. The default is native amino-
            acid k-mers only.

    Returns:
        Complete validated mapping for every observed feature key. Unsupported
        modalities conservatively contain positives only, making absence unknown.

    Raises:
        InputValidationError: If evidence references unknown proteins, explicit
            universes contradict positives, or universal type names are invalid.
    """

    known = _validated_protein_ids(values=protein_ids, context="campaign proteins")
    universal_types = _validated_feature_types(values=universally_assessed_feature_types)
    feature_proteins: dict[FeatureKey, set[str]] = defaultdict(set)
    universal_feature_keys: set[FeatureKey] = set()
    for feature in features:
        key = _validated_feature_key(value=(feature.feature_type, feature.feature_id))
        if feature.protein_id not in known:
            raise InputValidationError(
                f"Feature {key!r} references unknown protein {feature.protein_id!r}."
            )
        if feature.feature_type in universal_types:
            universal_feature_keys.add(key)
        else:
            feature_proteins[key].add(feature.protein_id)
    universes: dict[FeatureKey, set[str] | frozenset[str]] = {
        key: set(positives) for key, positives in feature_proteins.items()
    }
    proved_keys: set[FeatureKey] = set(universal_feature_keys)
    if isinstance(explicit_assessment_universes, Mapping):
        raise InputValidationError(
            "explicit_assessment_universes must be a collection of per-feature mappings."
        )
    for supplied_universes in explicit_assessment_universes:
        if not isinstance(supplied_universes, Mapping):
            raise InputValidationError(
                "Each explicit assessment-universe source must be a mapping."
            )
        explicit_positives = {key: feature_proteins.get(key, set()) for key in supplied_universes}
        explicit = normalise_feature_assessment_universes(
            universes=supplied_universes,
            known_protein_ids=known,
            feature_proteins=explicit_positives,
        )
        if explicit is None:  # pragma: no cover - guarded by the branch above
            raise InputValidationError("Explicit assessment universes were not normalised.")
        for key, assessed in explicit.items():
            if key in universal_feature_keys and assessed != known:
                raise InputValidationError(
                    f"Universally assessed feature {key!r} has an incomplete explicit "
                    "assessment universe."
                )
            if key in universes:
                universes[key].update(assessed)
                proved_keys.add(key)
    for key in universal_feature_keys:
        universes[key] = known
    assessed_domains: dict[str, set[str]] = defaultdict(set)
    for assessment in domain_assessments:
        if assessment.protein_id not in known:
            raise InputValidationError(
                f"Domain assessment references unknown protein {assessment.protein_id!r}."
            )
        if assessment.assessment_status in {
            DomainAssessmentStatus.ASSESSED_WITH_HIT,
            DomainAssessmentStatus.ASSESSED_NO_HIT,
        }:
            assessed_domains[assessment.domain_authority.casefold()].add(assessment.protein_id)
    for key in feature_proteins:
        feature_type, feature_id = key
        authority = feature_id.partition(":")[0].casefold()
        if (
            feature_type
            in {
                "PFAM_DOMAIN",
                "PFAM_ARCHITECTURE",
                "DOMAIN",
                "DOMAIN_ARCHITECTURE",
            }
            and authority in assessed_domains
        ):
            universes[key].update(assessed_domains[authority])
            proved_keys.add(key)
    conservative = sorted(set(universes) - proved_keys)
    if conservative:
        LOGGER.warning(
            "%d feature assessment universes contain observed positives only; "
            "unobserved values remain unknown until an assessment ledger is supplied.",
            len(conservative),
        )
    normalised = normalise_feature_assessment_universes(
        universes=universes,
        known_protein_ids=known,
        feature_proteins=feature_proteins,
    )
    if normalised is None:  # pragma: no cover - universes is always explicit
        raise InputValidationError("Feature assessment universes were not composed.")
    LOGGER.info(
        "Composed feature assessment universes definitions=%d "
        "universal_definitions=%d shared_universal_proteins=%d",
        len(normalised),
        len(universal_feature_keys),
        len(known),
    )
    return normalised


def normalise_feature_assessment_universes(
    *,
    universes: FeatureAssessmentUniverses | None,
    known_protein_ids: Collection[str],
    feature_proteins: Mapping[FeatureKey, Collection[str]],
) -> NormalisedAssessmentUniverses | None:
    """Validate and freeze per-feature assessment universes.

    A protein in a feature's universe was successfully assessed for that
    feature. A protein outside the universe is unknown, not an observed zero.
    Every observed positive must therefore be present in the matching universe.
    Repeated references to the same immutable universe are validated once and
    preserved by identity; positive collections are inspected without copying.

    Args:
        universes: Proteins successfully assessed for each feature, or ``None``
            only when every supplied feature is known to have been assessed for
            every relevant protein.
        known_protein_ids: Complete protein identifier universe for the stage.
        feature_proteins: Observed positive proteins for each feature.

    Returns:
        Frozen, validated universes, or ``None`` for the explicit legacy
        all-proteins-assessed contract.

    Raises:
        InputValidationError: If identifiers or feature keys are malformed, an
            assessment references an unknown protein, or an observed positive
            is not declared assessed.
    """

    known = _validated_protein_ids(values=known_protein_ids, context="known proteins")
    positives: dict[FeatureKey, Collection[str]] = {}
    for raw_key, raw_proteins in feature_proteins.items():
        key = _validated_feature_key(value=raw_key)
        _validate_protein_id_membership(
            values=raw_proteins,
            context=f"positive proteins for feature {key!r}",
            known_protein_ids=known,
        )
        positives[key] = raw_proteins
    if universes is None:
        return None
    normalised: NormalisedAssessmentUniverses = {}
    validated_collections: dict[int, tuple[Collection[str], frozenset[str]]] = {}
    for raw_key, raw_proteins in universes.items():
        key = _validated_feature_key(value=raw_key)
        if key in normalised:
            raise InputValidationError(f"Duplicate assessment universe for feature {key!r}.")
        cached = validated_collections.get(id(raw_proteins))
        if cached is not None and cached[0] is raw_proteins:
            protein_ids = cached[1]
        else:
            protein_ids = _validated_protein_ids(
                values=raw_proteins,
                context=f"assessment universe for feature {key!r}",
            )
            unknown = protein_ids - known
            if unknown:
                raise InputValidationError(
                    f"Assessment universe for feature {key!r} contains unknown proteins: "
                    f"{sorted(unknown)[:10]}"
                )
            validated_collections[id(raw_proteins)] = (raw_proteins, protein_ids)
        normalised[key] = protein_ids
    for key, positive_proteins in positives.items():
        if key not in normalised:
            raise InputValidationError(f"Missing assessment universe for observed feature {key!r}.")
        unassessed_positives = sorted(
            protein_id
            for protein_id in positive_proteins
            if protein_id not in normalised[key]
        )
        if unassessed_positives:
            raise InputValidationError(
                f"Observed positive proteins for feature {key!r} were not declared assessed: "
                f"{unassessed_positives[:10]}"
            )
    return normalised


def _validate_protein_id_membership(
    *,
    values: Collection[str],
    context: str,
    known_protein_ids: frozenset[str],
) -> None:
    """Validate protein identifiers without copying a large collection.

    Args:
        values: Candidate protein identifiers.
        context: Human-readable location used in validation errors.
        known_protein_ids: Complete permitted protein universe.

    Raises:
        InputValidationError: If the collection is malformed or contains an
            identifier outside the known universe.
    """

    if isinstance(values, (str, bytes)) or not isinstance(values, Collection):
        raise InputValidationError(f"{context.capitalize()} must be a collection of strings.")
    unknown: list[str] = []
    for value in values:
        if not isinstance(value, str) or not value.strip():
            raise InputValidationError(
                f"{context.capitalize()} contains empty or non-string protein identifiers."
            )
        if value not in known_protein_ids and len(unknown) < 10:
            unknown.append(value)
    if unknown:
        raise InputValidationError(
            f"Feature {context.removeprefix('positive proteins for feature ')} contains "
            f"proteins outside the known universe: {sorted(unknown)}"
        )


def _validated_feature_key(*, value: object) -> FeatureKey:
    """Return one canonical feature key after strict validation.

    Args:
        value: Candidate ``(feature_type, feature_id)`` pair.

    Returns:
        Validated feature key.

    Raises:
        InputValidationError: If the key is not two non-empty strings.
    """

    if (
        not isinstance(value, tuple)
        or len(value) != 2
        or not all(isinstance(item, str) and item.strip() for item in value)
    ):
        raise InputValidationError(
            "Feature-assessment keys must be (feature_type, feature_id) tuples "
            "of non-empty strings."
        )
    return value


def _validated_feature_types(*, values: Collection[str]) -> frozenset[str]:
    """Return validated universally assessed feature-type names.

    Args:
        values: Candidate feature-type collection.

    Returns:
        Frozen feature-type names.

    Raises:
        InputValidationError: If the input is scalar or contains invalid names.
    """

    if isinstance(values, (str, bytes)) or not isinstance(values, Collection):
        raise InputValidationError("Universal feature types must be a collection of strings.")
    output = frozenset(values)
    if any(not isinstance(value, str) or not value.strip() for value in output):
        raise InputValidationError("Universal feature types must contain non-empty strings.")
    return output


def _validated_protein_ids(*, values: Collection[str], context: str) -> frozenset[str]:
    """Return validated protein identifiers as an immutable set.

    Args:
        values: Candidate collection of identifiers.
        context: Human-readable location used in validation errors.

    Returns:
        Frozen protein identifier set.

    Raises:
        InputValidationError: If the value is a scalar string or contains an
            empty or non-string identifier.
    """

    if isinstance(values, (str, bytes)) or not isinstance(values, Collection):
        raise InputValidationError(f"{context.capitalize()} must be a collection of strings.")
    output = frozenset(values)
    invalid = [value for value in output if not isinstance(value, str) or not value.strip()]
    if invalid:
        raise InputValidationError(
            f"{context.capitalize()} contains empty or non-string protein identifiers."
        )
    return output
