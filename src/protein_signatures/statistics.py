"""Dependency-light exact statistics for categorical signature features."""

from __future__ import annotations

import math

from .errors import InputValidationError


def fisher_exact_two_sided(*, a: int, b: int, c: int, d: int) -> tuple[float | None, float]:
    """Calculate an odds ratio and two-sided Fisher exact p-value.

    Args:
        a: Target inference units with the feature.
        b: Target inference units without the feature.
        c: Background inference units with the feature.
        d: Background inference units without the feature.

    Returns:
        Odds ratio (or ``None`` when undefined) and exact p-value.

    Raises:
        InputValidationError: If counts are negative or non-integral.
    """

    counts = (a, b, c, d)
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 0 for value in counts):
        raise InputValidationError("Fisher exact counts must be non-negative integers.")
    numerator = a * d
    denominator = b * c
    if denominator == 0:
        odds_ratio = math.inf if numerator > 0 else None
    else:
        odds_ratio = numerator / denominator
    total = sum(counts)
    if total == 0:
        return odds_ratio, 1.0
    row_target = a + b
    column_present = a + c
    minimum = max(0, row_target - (total - column_present))
    maximum = min(row_target, column_present)
    observed = _hypergeometric_probability(
        x=a,
        row_target=row_target,
        column_present=column_present,
        total=total,
    )
    p_value = sum(
        probability
        for x in range(minimum, maximum + 1)
        if (
            probability := _hypergeometric_probability(
                x=x,
                row_target=row_target,
                column_present=column_present,
                total=total,
            )
        )
        <= observed * (1.0 + 1e-12)
    )
    return odds_ratio, min(1.0, p_value)


def benjamini_hochberg(*, p_values: tuple[float, ...]) -> tuple[float, ...]:
    """Adjust p-values using the Benjamini-Hochberg FDR procedure.

    Args:
        p_values: Finite p-values from zero to one.

    Returns:
        Adjusted q-values in original order.

    Raises:
        InputValidationError: If a p-value is invalid.
    """

    if any(not math.isfinite(value) or not 0.0 <= value <= 1.0 for value in p_values):
        raise InputValidationError("P-values must be finite values from 0.0 to 1.0.")
    count = len(p_values)
    if count == 0:
        return ()
    ordered = sorted(enumerate(p_values), key=lambda item: (item[1], item[0]))
    adjusted = [1.0] * count
    running = 1.0
    for reverse_rank, (index, p_value) in enumerate(reversed(ordered), start=1):
        rank = count - reverse_rank + 1
        running = min(running, p_value * count / rank)
        adjusted[index] = min(1.0, running)
    return tuple(adjusted)


def wilson_score_interval(
    *, successes: int, trials: int, z_score: float = 1.959963984540054
) -> tuple[float, float]:
    """Calculate a Wilson score interval for a binomial proportion.

    Args:
        successes: Number of units carrying the feature.
        trials: Number of assessed independent units.
        z_score: Positive normal quantile. The default gives a two-sided 95% interval.

    Returns:
        Lower and upper bounds constrained to zero through one.

    Raises:
        InputValidationError: If counts or the quantile are invalid.
    """

    if (
        not isinstance(successes, int)
        or isinstance(successes, bool)
        or not isinstance(trials, int)
        or isinstance(trials, bool)
        or trials < 1
        or successes < 0
        or successes > trials
        or not math.isfinite(z_score)
        or z_score <= 0.0
    ):
        raise InputValidationError(
            "Wilson interval requires 0 <= successes <= trials, trials >= 1, "
            "and a positive finite z-score."
        )
    proportion = successes / trials
    squared = z_score * z_score
    denominator = 1.0 + squared / trials
    centre = (proportion + squared / (2.0 * trials)) / denominator
    half_width = (
        z_score
        * math.sqrt(proportion * (1.0 - proportion) / trials + squared / (4.0 * trials * trials))
        / denominator
    )
    return max(0.0, centre - half_width), min(1.0, centre + half_width)


def _hypergeometric_probability(
    *, x: int, row_target: int, column_present: int, total: int
) -> float:
    """Calculate one fixed-margin hypergeometric table probability.

    Args:
        x: Target-with-feature count.
        row_target: Target row total.
        column_present: Feature-present column total.
        total: Overall sample size.

    Returns:
        Probability of the table.
    """

    if x < 0 or x > row_target or x > column_present:
        return 0.0
    other = row_target - x
    absent_total = total - column_present
    if other < 0 or other > absent_total:
        return 0.0
    log_probability = (
        _log_combination(n=column_present, k=x)
        + _log_combination(n=absent_total, k=other)
        - _log_combination(n=total, k=row_target)
    )
    return math.exp(log_probability)


def _log_combination(*, n: int, k: int) -> float:
    """Return the log of an exact integer combination.

    Args:
        n: Population size.
        k: Selected count.

    Returns:
        Natural logarithm of ``n choose k``.
    """

    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
