"""
Client-side QC aggregation functions.

Clients aggregate QC data locally to compute exclusion lists without the server ever
seeing plaintext statistics.
"""

from __future__ import annotations

import logging
from typing import Set

import numpy as np
from scipy.stats import chi2

# Note: chi2.sf(x, df) = 1 - chi2.cdf(x, df) but more numerically stable for small p-values

# Module-level logger for QC aggregation
logger = logging.getLogger(__name__)


def _compute_exclusion_list(
    counts_sum: np.ndarray,
    missing_sum: np.ndarray,
    thresholds: np.ndarray,
    debug: bool = False,
) -> Set[int]:
    """
    Compute exclusion list from aggregated QC statistics.

    Args:
        counts_sum: Aggregated genotype counts, shape (n_snps, 3)
        missing_sum: Aggregated missing data, shape (n_snps, 2)
        thresholds: Threshold values [maf, missing, hwe]
        debug: If True, log debug information

    Returns:
        Set of SNP indices to exclude
    """
    if counts_sum.size == 0 or missing_sum.size == 0:
        if debug:
            logger.debug(
                f"[QC DEBUG] Empty arrays: counts_sum.size={counts_sum.size}, missing_sum.size={missing_sum.size}"
            )
        return set()

    maf_final = float(thresholds[0]) if len(thresholds) > 0 else 0.01
    miss_final = float(thresholds[1]) if len(thresholds) > 1 else 0.1
    hwe_final = float(thresholds[2]) if len(thresholds) > 2 else 1e-6

    if debug:
        logger.debug(
            f"[QC DEBUG] Thresholds: MAF={maf_final}, Missing={miss_final}, HWE={hwe_final}"
        )

    n_snps = counts_sum.shape[0]
    exclude_indices: Set[int] = set()

    excluded_by_maf = 0
    excluded_by_missing = 0
    excluded_by_hwe = 0
    excluded_by_zero_geno = 0

    for i in range(n_snps):
        # MAF
        N_AA, N_Aa, N_aa = counts_sum[i]
        N_obs, N_miss = missing_sum[i]

        total_geno = N_AA + N_Aa + N_aa
        if total_geno == 0:
            exclude_indices.add(i)
            excluded_by_zero_geno += 1
            if debug and excluded_by_zero_geno <= 3:
                logger.debug(f"[QC DEBUG] SNP {i}: Excluded by zero genotypes")
            continue

        freqA = (2 * N_AA + N_Aa) / (2 * total_geno)
        maf = min(freqA, 1.0 - freqA)
        if maf < maf_final:
            exclude_indices.add(i)
            excluded_by_maf += 1
            if debug and excluded_by_maf <= 3:
                logger.debug(
                    f"[QC DEBUG] SNP {i}: Excluded by MAF (maf={maf:.6f} < {maf_final})"
                )
            continue

        # Missing rate (per-SNP)
        denom = N_obs + N_miss
        if denom == 0:
            exclude_indices.add(i)
            excluded_by_zero_geno += 1
            if debug and excluded_by_zero_geno <= 3:
                logger.debug(f"[QC DEBUG] SNP {i}: Excluded by zero denominator")
            continue
        missing_rate = N_miss / denom
        if missing_rate > miss_final:
            exclude_indices.add(i)
            excluded_by_missing += 1
            if debug and excluded_by_missing <= 3:
                logger.debug(
                    f"[QC DEBUG] SNP {i}: Excluded by missing rate ({missing_rate:.6f} > {miss_final})"
                )
            continue

        # chi-squared based HWE test
        # Only perform HWE test if we have enough observations (at least 5 genotypes)
        if total_geno < 5:
            # Skip HWE test for SNPs with too few observations
            if debug and excluded_by_hwe <= 3:
                logger.debug(
                    f"[QC DEBUG] SNP {i}: Skipping HWE test (too few observations: {total_geno})"
                )
        else:
            p = freqA
            q = 1.0 - p
            # Avoid edge cases where p or q is exactly 0 or 1
            if p <= 0 or p >= 1 or q <= 0 or q >= 1:
                # Skip HWE test for monomorphic SNPs
                if debug and excluded_by_hwe <= 3:
                    logger.debug(
                        f"[QC DEBUG] SNP {i}: Skipping HWE test (monomorphic: p={p:.6f})"
                    )
            else:
                E_AA = p * p * total_geno
                E_Aa = 2.0 * p * q * total_geno
                E_aa = q * q * total_geno

                # Only compute chi-squared if all expected counts are reasonable
                if E_AA < 0.5 or E_Aa < 0.5 or E_aa < 0.5:
                    # Skip HWE test if expected counts are too small
                    if debug and excluded_by_hwe <= 3:
                        logger.debug(
                            f"[QC DEBUG] SNP {i}: Skipping HWE test (expected counts too small: E_AA={E_AA:.2f}, E_Aa={E_Aa:.2f}, E_aa={E_aa:.2f})"
                        )
                else:
                    # Special case: When N_aa = 0 but E_aa is large, this indicates severe HWE violation
                    # However, PLINK may handle this case differently (e.g., skip HWE test or use exact test)
                    # For now, we'll still perform the chi-squared test but note this pattern
                    if N_aa == 0 and E_aa > 5.0:
                        # This is a severe HWE violation: expected homozygous recessive but observed none
                        # PLINK typically excludes such SNPs, so we do the same
                        if debug and excluded_by_hwe <= 3:
                            logger.debug(
                                f"[QC DEBUG] SNP {i}: Severe HWE violation (N_aa=0 but E_aa={E_aa:.2f})"
                            )

                    chi_sq = 0.0
                    for obs, exp in [(N_AA, E_AA), (N_Aa, E_Aa), (N_aa, E_aa)]:
                        if exp > 1e-9:
                            chi_sq += (obs - exp) ** 2 / exp

                    # Compute p-value, handling edge cases
                    if chi_sq > 0:
                        # Use survival function (1 - CDF) for better numerical stability
                        # For very large chi_sq, pval will be very small but not exactly 0
                        pval_hwe = chi2.sf(
                            chi_sq, df=1
                        )  # sf = survival function = 1 - cdf
                        # Clamp p-value to avoid numerical issues
                        pval_hwe = max(1e-300, min(1.0, pval_hwe))
                    else:
                        # Perfect HWE, p-value = 1.0
                        pval_hwe = 1.0

                    if debug and excluded_by_hwe <= 3:
                        logger.debug(
                            f"[QC DEBUG] SNP {i}: HWE test: chi_sq={chi_sq:.4f}, pval={pval_hwe:.2e}, obs=[{N_AA},{N_Aa},{N_aa}], exp=[{E_AA:.2f},{E_Aa:.2f},{E_aa:.2f}]"
                        )

                    if pval_hwe < hwe_final:
                        exclude_indices.add(i)
                        excluded_by_hwe += 1
                        if debug and excluded_by_hwe <= 3:
                            logger.debug(
                                f"[QC DEBUG] SNP {i}: Excluded by HWE (p={pval_hwe:.2e} < {hwe_final})"
                            )

    if debug:
        logger.debug(
            f"[QC DEBUG] Exclusion summary: MAF={excluded_by_maf}, Missing={excluded_by_missing}, HWE={excluded_by_hwe}, ZeroGeno={excluded_by_zero_geno}, Total={len(exclude_indices)}"
        )

    return exclude_indices
