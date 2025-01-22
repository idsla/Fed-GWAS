## Architecture

<!-- include the architecture charts -->
The Federated GWAS Pipeline is designed to enable collaborative genomic data analysis across multiple institutions without sharing raw data, thereby preserving privacy and complying with data protection regulations. The architecture is built around the principles of decentralized processing, secure communication, and centralized aggregation of results.

![Alt Text](/images/Architecture.png)

## Pipeline Flow

<!-- include the flow charts -->
This flow diagram outlines the process of federated learning in the GWAS pipeline:

1. Local data is cleaned and processed at each client.
2. Summary statistics (like kinship coefficients) are securely shared with the server.
3. The server aggregates the statistics to create global models, which are then sent back to the clients for further local analysis.

# Components

## 1. Quality Control

## Components
**Quality Control (QC)** module ensures that genomic and phenotype data are clean and reliable before further analysis. Additionally, it computes kinship coefficients to measure genetic relatedness between individuals. The following steps are performed:

1. Quality Control
- **Missing Rate Calculation**: Identifies missing values in genotype data to determine call rates for each SNP.
- **Minor Allele Frequency (MAF) Calculation**: Ensures that SNPs with rare alleles are identified and optionally filtered.
- **Hardy-Weinberg Equilibrium Test**: Checks for deviations from Hardy-Weinberg equilibrium to remove unreliable SNPs.
- **Sex Check**: Ensures the genetic sex of each individual aligns with reported phenotypes.
- **Filtering**: Applies threshold-based filters (e.g., MAF thresholds, call rate thresholds) to clean the dataset.
- **Kinship Matrix Calculation**: Computes kinship coefficients locally using SNP data. These coefficients indicate the genetic relationship between individuals.

<!-- add introduction to quality control module including kinship -->
## 2. Association Analysis

2. Association Analysis
**Association Analysis** module examines the association between SNPs and phenotypes across the population. Each client performs the following steps:

<!-- add introduction to association analysis module -->
- **Logistic Regression**: Performs logistic regression for binary phenotypes (e.g., disease status).
- **Linear Regression**: Conducts linear regression for continuous phenotypes (e.g., height, weight).
- **Chi-Square Test**: Tests for significant associations between individual SNPs and phenotypes using a chi-square test.
- **Global P-Value Aggregation**: Clients share p-values with the server, which aggregates them to generate global association results.