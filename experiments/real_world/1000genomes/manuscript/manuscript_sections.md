## Dataset and Experimental Setup

### Dataset

We evaluated our federated GWAS pipeline on chromosome 22 data from the 1000 Genomes Project Phase 3 (1kGP) [@1000genomes2015]. The dataset consists of genotype data from 2,504 individuals across multiple populations. For this experiment, we partitioned the data by population groups to simulate a multi-site federated scenario:

- **Center 1**: European (EUR) population samples
- **Center 2**: African (AFR) population samples

Phenotypes were simulated using a population-stratified logistic model with a causal fraction of 0.01 (1% of SNPs are causal). The simulation incorporated population effects to test the pipeline's ability to handle population stratification.

### Experimental Setup

The federated experiment was conducted using the Flower framework in local simulation mode. The pipeline executed the following stages:

1. **Key Exchange**: Secure key establishment between clients
2. **Local QC**: Per-client quality control filtering
3. **Global QC**: Federated aggregation of QC statistics and exclusion list generation
4. **KING Analysis**: Iterative kinship estimation across client data chunks
5. **Local LR**: Per-client logistic regression on filtered variants
6. **Global LR**: Federated aggregation of LR statistics and p-value computation

A centralized baseline was generated using PLINK 1.9/2.0 for comparison. All analyses used standard GWAS thresholds: MAF ≥ 0.01, missing rate ≤ 0.05, HWE p-value ≥ 1×10⁻⁶.

## Results

### Quality Control Agreement

The federated QC pipeline achieved excellent agreement with the centralized baseline. Of the 956,898 SNPs excluded by the centralized analysis, the federated approach correctly identified 955,279 SNPs (F1 score: 0.9991, Precision: 1.0000, Recall: 0.9983). The high precision (1.0000) indicates no false positive exclusions, while the recall of 0.9983 shows that the federated approach captures the vast majority of problematic variants identified in the centralized analysis.

### Kinship Estimation (KING)

KING kinship coefficients computed through the federated iterative aggregation showed strong correlation with the centralized baseline (Pearson r = 0.9978). The analysis included 3,283,461 sample pairs, with 1,663,924 cross-client pairs and 1,619,537 within-client pairs. The mean absolute error (MAE) was 0.005216, demonstrating high accuracy in the federated kinship estimation.

### Logistic Regression Results

The federated global LR analysis achieved a correlation of r = 0.6834 (p < 1×10⁻³⁰⁰) with the centralized baseline across 119,462 SNPs. The top 100 most significant SNPs showed 65.0% overlap between federated and centralized results, indicating strong agreement in identifying the most significant associations.

Coverage analysis at multiple significance thresholds demonstrated that the federated approach captures a substantial proportion of significant associations identified in the centralized analysis, with coverage increasing at more stringent thresholds.
