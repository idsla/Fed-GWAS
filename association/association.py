import pandas as pd
import numpy as np
import statsmodels.api as sm

def linear_regression_snp(snp, phenotype_):
    if snp.isnull().any() or phenotype_.isnull().any():
        raise ValueError("Missing values detected in SNP or phenotype data.")
    
    if snp.var() == 0:
        raise ValueError("SNP has zero variance.")
    
    x = sm.add_constant(snp)
    model = sm.OLS(phenotype_, x)
    result = model.fit()
    return result.pvalues.iloc[1]

def snpwise_linear_regression(genotype_df, phenotype_):
    p_values = []
    for snp in genotype_df.columns:
        try:
            p_value = linear_regression_snp(genotype_df[snp], phenotype_)
            p_values.append(p_value)
        except ValueError as ve:
            if "zero variance" in str(ve):
                print(f"Skipping SNP {snp} due to zero variance.")
                p_values.append(np.nan)
            else:
                print(f"Error processing SNP {snp}: {ve}")
                p_values.append(np.nan)
        except Exception as e:
            print(f"Unexpected error processing SNP {snp}: {e}")
            p_values.append(np.nan)
    return pd.Series(p_values, index=genotype_df.columns)