# Correctness Experiment: Tiny Scale, Even Partition

## Purpose

This experiment validates that federated GWAS results match centralized PLINK results. It serves as a baseline correctness test for the federated pipeline.

## Configuration

- **Scale**: Tiny (500 samples, 5,000 SNPs)
- **Clients**: 2
- **Partition Strategy**: Even (equal samples per client)
- **Expected Runtime**: 5-10 minutes

## Setup

### 1. Generate Data

```bash
python pipeline/simulation/simulated_data/generate_synthetic_data.py \
  --scale tiny \
  --partition-strategy even \
  --seed 42 \
  --output-dir experiments/correctness/tiny_even/data
```

### 2. Generate Baseline

```bash
python experiments/tools/generate_baseline.py \
  experiments/correctness/tiny_even/config.yaml
```

This will create a centralized PLINK baseline in `data/tiny/centralized_baseline/` for comparison.

### 3. Run Experiment

Run the experiment via the Flower CLI using `experiments/correctness/tiny_even/config.yaml`.

## Analysis

After the experiment completes:

### 1. Collect Metrics

```bash
python experiments/tools/evaluation/metrics_collector.py \
  experiments/correctness/tiny_even/results
```

### 2. Compare with Baseline

```bash
# QC agreement
python experiments/tools/evaluation/qc/qc_evaluator.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline

# LR agreement
python experiments/tools/evaluation/lr/lr_evaluator.py \
  experiments/correctness/tiny_even/results \
  --baseline experiments/correctness/tiny_even/data/tiny/centralized_baseline
```

This will generate:
- `qc_report.md` - QC agreement report
- `lr_report.md` - LR agreement report
- `manhattan_plot.png` - Manhattan plot comparison
- `lr_correlation_plots.png` - LR correlation plots

## Expected Results

- **QC Agreement**: F1 score > 0.95
- **LR Correlation**: Pearson r > 0.99
- **Top 100 SNPs Overlap**: > 80%

## Success Criteria

The experiment passes if:
- All stages complete successfully
- QC F1 score > 0.95
- LR correlation > 0.99

## Files Generated

- `results/experiment_report.md` - Experiment summary
- `results/collected_metrics.json` - Aggregated metrics
- `results/qc_report.md` - QC comparison with baseline
- `results/lr_report.md` - LR comparison with baseline
- `results/manhattan_plot.png` - Visualization
- `results/lr_correlation_plots.png` - LR correlation visualizations
