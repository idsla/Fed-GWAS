import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from scipy import stats
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import warnings
from pysnptools.snpreader import Bed
import argparse

# File paths
parser = argparse.ArgumentParser(description="Perform QC on GWAS data.")
parser.add_argument('--fam', required=True, help="Path to the .fam file")
parser.add_argument('--bim', required=True, help="Path to the .bim file")
parser.add_argument('--bed', required=True, help="Path to the .bed file")
parser.add_argument('--threshold', type=float, default=0.05)
args = parser.parse_args()

bed_path = args.bed
bim_path = args.bim
fam_path = args.fam

# Read genotype data using pysnptools
snp_reader = Bed(bed_path, count_A1=False)
geno = snp_reader.read()

# Extract genotype data (X) and phenotype data (y)
X = geno.val  # Genotype data
y = np.loadtxt(fam_path, usecols=[5])  # Assuming phenotype is in the 6th column of .fam file

# Convert PHENOTYPE values: 2 -> 1 (case), 1 -> 0 (control)
y = np.where(y == 2, 1, 0)

# Check the distribution of the phenotype data
unique_classes, counts = np.unique(y, return_counts=True)
print("Class distribution in phenotype data:", dict(zip(unique_classes, counts)))

# If there is only one class, raise an exception
if len(unique_classes) < 2:
    raise ValueError("The phenotype data contains only one class. Logistic regression requires at least two classes.")

# Read SNP IDs from the .bim file
snp_ids = pd.read_csv(bim_path, sep='\s+', header=None)[1].values

# Handle missing values using mean imputation
imputer = SimpleImputer(strategy='mean')
X_imputed = imputer.fit_transform(X)

# Normalize data
scaler = StandardScaler()
X_normalized = scaler.fit_transform(X_imputed)

# Function to perform logistic regression using sklearn and statsmodels
def logistic_regression_comparison(X, y):
    p_values_sklearn = []
    p_values_statsmodels = []
    num_snps = X.shape[1]

    for snp_idx in range(num_snps):
        snp_data = X[:, snp_idx].reshape(-1, 1)

        # Sklearn logistic regression
        try:
            model_sklearn = LogisticRegression(solver='liblinear')
            model_sklearn.fit(snp_data, y)
            
            # Extract the model parameters
            coef = model_sklearn.coef_[0][0]
            intercept = model_sklearn.intercept_[0]

            # Manually compute standard errors and p-values
            pred_probs = model_sklearn.predict_proba(snp_data)[:, 1]
            X_design = np.hstack([np.ones((snp_data.shape[0], 1)), snp_data])
            V = np.diagflat(pred_probs * (1 - pred_probs))
            cov_matrix = np.linalg.inv(X_design.T @ V @ X_design)
            se = np.sqrt(np.diag(cov_matrix))
            z_scores = coef / se[1]
            p_value_sklearn = 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
            p_values_sklearn.append(p_value_sklearn)
        except Exception as e:
            p_values_sklearn.append(np.nan)

        # Statsmodels logistic regression
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=ConvergenceWarning)
                model_statsmodels = sm.Logit(y, sm.add_constant(snp_data)).fit(disp=False)
            p_value_statsmodels = model_statsmodels.pvalues[1]  # Get p-value of the SNP coefficient
            p_values_statsmodels.append(p_value_statsmodels)
        except Exception as e:
            p_values_statsmodels.append(np.nan)         

    return np.array(p_values_sklearn), np.array(p_values_statsmodels)

# Run logistic regression comparison
p_values_sklearn_log, p_values_statsmodels_log = logistic_regression_comparison(X_normalized, y)

# Remove SNPs with NaN p-values to ensure lengths match
""" valid_indices = ~np.isnan(p_values_sklearn_log) & ~np.isnan(p_values_statsmodels_log)
snp_ids = snp_ids[valid_indices]
p_values_sklearn_log = p_values_sklearn_log[valid_indices]
p_values_statsmodels_log = p_values_statsmodels_log[valid_indices] """

# Function to perform linear regression using sklearn and statsmodels
def linear_regression_comparison(X, y):
    p_values_sklearn = []
    p_values_statsmodels = []
    num_snps = X.shape[1]

    for snp_idx in range(num_snps):
        snp_data = X[:, snp_idx].reshape(-1, 1)

        # Sklearn linear regression
        try:
            model_sklearn = LinearRegression()
            model_sklearn.fit(snp_data, y)
            coef = model_sklearn.coef_[0]
            intercept = model_sklearn.intercept_
            y_pred = model_sklearn.predict(snp_data)
            residuals = y - y_pred
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            r2 = 1 - (ss_res / ss_tot)
            var_b = np.var(residuals) / (np.var(snp_data) * len(y))
            se_b = np.sqrt(var_b)
            if se_b == 0:
                p_values_sklearn.append(np.nan)
            else:
                t_stat = coef / se_b
                p_value_sklearn = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=len(y) - 2))
                p_values_sklearn.append(p_value_sklearn)
        except Exception as e:
            p_values_sklearn.append(np.nan)

        # Statsmodels linear regression
        try:
            model_statsmodels = sm.OLS(y, sm.add_constant(snp_data)).fit()
            p_value_statsmodels = model_statsmodels.pvalues[1]  # Get p-value of the SNP coefficient
            p_values_statsmodels.append(p_value_statsmodels)
        except Exception as e:
            p_values_statsmodels.append(np.nan)

    return np.array(p_values_sklearn), np.array(p_values_statsmodels)

# Run linear regression comparison
p_values_sklearn_lin, p_values_statsmodels_lin = linear_regression_comparison(X_normalized, y)

# Compare p-values
comparison_df = pd.DataFrame({
    'SNP ID': snp_ids,
    'p-value Statsmodels Logistic': p_values_statsmodels_log,
    'p-value Sklearn Logistic': p_values_sklearn_log,
    'p-value Statsmodels Linear': p_values_statsmodels_lin,
    'p-value Sklearn Linear': p_values_sklearn_lin,
    
})

# Set the significance level
threshold = args.threshold
num_tests = len(snp_ids)  # Total number of SNPs

# Calculate Bonferroni-corrected significance threshold
alpha_bonferroni = threshold/num_tests

# Apply Bonferroni correction and determine significance
comparison_df['Significant Statsmodels Logistic'] = comparison_df['p-value Statsmodels Logistic'] < alpha_bonferroni
comparison_df['Significant Sklearn Logistic'] = comparison_df['p-value Sklearn Logistic'] < alpha_bonferroni
comparison_df['Significant Statsmodels Linear'] = comparison_df['p-value Statsmodels Linear'] < alpha_bonferroni
comparison_df['Significant Sklearn Linear'] = comparison_df['p-value Sklearn Linear'] < alpha_bonferroni

# Save comparison results to a CSV file with significance information
comparison_df.to_csv('comparison_regression_results_with_significance.csv', index=False)
# Print the number of significant SNPs for each method
num_significant_statsmodels_log = comparison_df['Significant Statsmodels Logistic'].sum()
num_significant_sklearn_log = comparison_df['Significant Sklearn Logistic'].sum()
num_significant_statsmodels_lin = comparison_df['Significant Statsmodels Linear'].sum()
num_significant_sklearn_lin = comparison_df['Significant Sklearn Linear'].sum()

print(f"Number of significant SNPs (Statsmodels Logistic): {num_significant_statsmodels_log}")
print(f"Number of significant SNPs (Sklearn Logistic): {num_significant_sklearn_log}")
print(f"Number of significant SNPs (Statsmodels Linear): {num_significant_statsmodels_lin}")
print(f"Number of significant SNPs (Sklearn Linear): {num_significant_sklearn_lin}")
print(f"Comparison results with significance saved to 'comparison_regression_results_with_significance.csv'.")
