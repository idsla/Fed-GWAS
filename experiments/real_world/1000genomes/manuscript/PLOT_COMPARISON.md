# Plot Comparison: Previous vs New Manuscript Plots

## Previous Plot (from evaluation tool)

**File**: `experiments/real_world/1000genomes/results/manhattan_plot.png`
- **Generated**: February 22, 2026
- **Size**: 455 KB
- **Resolution**: 1784 × 1181 pixels (150 DPI)
- **Format**: Side-by-side comparison (Federated top, Baseline bottom)
- **Data Source**: 
  - Federated: `lr_results_latest.assoc.logistic` (120,806 SNPs)
  - Baseline: `lr.assoc.logistic` (144,369 SNPs)
- **Purpose**: Quick comparison plot from evaluation tool

## New Manuscript Plots

### Individual Plots (300 DPI, publication-ready)

1. **`manhattan_federated.png`** (696 KB, 3570 × 1770 pixels)
   - Federated GWAS results only
   - High-resolution for publication
   - Includes genome-wide significance line

2. **`manhattan_baseline.png`** (773 KB, 3570 × 1770 pixels)
   - Centralized baseline results only
   - High-resolution for publication
   - Includes genome-wide significance line

3. **`qq_federated.png`** (161 KB)
   - Q-Q plot for federated results
   - Shows p-value distribution and genomic inflation factor (λ)

4. **`qq_baseline.png`** (160 KB)
   - Q-Q plot for baseline results
   - Shows p-value distribution and genomic inflation factor (λ)

### Comparison Plot (300 DPI, publication-ready)

5. **`manhattan_comparison.png`** (1.5 MB, high resolution)
   - Side-by-side comparison (Federated top, Baseline bottom)
   - **Matches the format of the previous plot** but at publication quality
   - Aligned by genomic position (CHR + BP) since SNP IDs differ
   - 119,465 SNPs matched between federated and baseline
   - Includes genome-wide significance lines on both panels

## Key Differences

### Data Matching
- **Previous plot**: Likely used de-anonymized SNP IDs or position matching
- **New comparison plot**: Uses CHR + BP position matching (119,465 SNPs matched)
- **Reason**: Federated results use anonymized SNP IDs (e.g., `rs4133873035`) while baseline uses position-based IDs (e.g., `22:16050654`)

### Resolution
- **Previous**: 150 DPI (adequate for screen viewing)
- **New**: 300 DPI (publication quality)

### Format Options
- **Previous**: Only side-by-side comparison
- **New**: Both individual plots (for separate figures) and comparison plot (for direct comparison)

## Recommendation for Manuscript

1. **For main results figure**: Use `manhattan_comparison.png` (side-by-side format, matches previous style)
2. **For supplementary figures**: Use individual plots (`manhattan_federated.png`, `manhattan_baseline.png`)
3. **For Q-Q plots**: Use both `qq_federated.png` and `qq_baseline.png` to show p-value distributions

## Data Validation

- ✅ All plots generated successfully
- ✅ Position-based matching works correctly (119,465 SNPs matched)
- ✅ P-value ranges are consistent between plots
- ✅ Genome-wide significance thresholds properly displayed
- ✅ High-resolution output suitable for publication
