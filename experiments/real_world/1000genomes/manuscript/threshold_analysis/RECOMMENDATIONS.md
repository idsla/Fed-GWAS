# Federated Screening Threshold Optimization

## Current Configuration
- **Current threshold**: p = 0.3
- **Coverage**: 100% (100/100 baseline significant SNPs preserved)
- **Federated significant SNPs**: 102,417
- **False positives**: 102,317

## Key Finding

**All tested thresholds (0.05 - 1.0) achieve 100% coverage** for genome-wide significant SNPs (p < 5×10⁻⁸). This means we can optimize for **filtering efficiency** (reducing false positives) without losing any true positives.

## Optimal Threshold Recommendation

### Best Option: **p = 0.1** (or even **p = 0.05**)

**Threshold p = 0.1:**
- ✅ **Coverage**: 100% (100/100 preserved)
- ✅ **Federated significant**: 57,946 SNPs
- ✅ **False positives**: 57,846 SNPs
- ✅ **Filtering ratio**: 0.0017
- **Improvement**: ~44% reduction in false positives compared to p = 0.3

**Threshold p = 0.05:**
- ✅ **Coverage**: 100% (100/100 preserved)
- ✅ **Federated significant**: 35,501 SNPs
- ✅ **False positives**: 35,401 SNPs
- ✅ **Filtering ratio**: 0.0028
- **Improvement**: ~65% reduction in false positives compared to p = 0.3

## Comparison Table

| Threshold | Coverage | Preserved | Federated Significant | False Positives | Filtering Ratio |
|-----------|----------|-----------|----------------------|-----------------|-----------------|
| **0.05**  | **1.0000** | **100** | **35,501** | **35,401** | **0.0028** |
| **0.1**   | **1.0000** | **100** | **57,946** | **57,846** | **0.0017** |
| 0.15      | 1.0000 | 100 | 74,497 | 74,397 | 0.0013 |
| 0.2       | 1.0000 | 100 | 86,930 | 86,830 | 0.0012 |
| 0.25      | 1.0000 | 100 | 95,428 | 95,328 | 0.0010 |
| **0.3** (current) | **1.0000** | **100** | **102,417** | **102,317** | **0.0010** |
| 0.35      | 1.0000 | 100 | 108,233 | 108,133 | 0.0009 |
| 0.4       | 1.0000 | 100 | 112,128 | 112,028 | 0.0009 |

## Recommendation

**Decrease the threshold from 0.3 to 0.1** (or even 0.05) to achieve:
- ✅ **100% coverage** (no missing true positives)
- ✅ **Better filtering** (fewer false positives)
- ✅ **More efficient downstream analysis** (fewer SNPs to process)

The threshold of **p = 0.1** provides an excellent balance:
- Maintains perfect coverage
- Reduces false positives by ~44% compared to current
- Still captures enough SNPs for robust federated analysis

For even more aggressive filtering, **p = 0.05** reduces false positives by ~65% while still maintaining 100% coverage.

## Implementation

To use the optimized threshold, update the federated screening threshold parameter in the pipeline configuration. The current threshold of 0.3 can be safely reduced to 0.1 (or 0.05) without any loss of coverage for genome-wide significant SNPs.
