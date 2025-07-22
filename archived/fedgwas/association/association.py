import numpy as np
import pandas as pd
import warnings
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from scipy import stats
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from pysnptools.snpreader import Bed
from fedgwas.io.reader import IODataHandler
import logging
logging.basicConfig(filename='report_association.log', filemode='w', level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')
class GWASAssociation:
    """
    A class to perform quality control and statistical analysis on GWAS data using both
    logistic and linear regression models from sklearn and statsmodels.

    Attributes:
    -----------
    threshold : float
        Significance level threshold for statistical tests.

    Methods:
    --------
    logistic_regression_comparison() -> pd.DataFrame:
        Performs logistic regression using sklearn and statsmodels and compares p-values.

    linear_regression_comparison() -> pd.DataFrame:
        Performs linear regression using sklearn and statsmodels and compares p-values.

    save_results(logistic_df: pd.DataFrame, linear_df: pd.DataFrame) -> None:
        Saves the results of logistic and linear regression comparisons to a CSV file.
    """

    def __init__(self, bed_path: str, bim_path: str, fam_path: str, threshold: float):
        """

        Initializes the GWASQC class with file paths for genotype and phenotype data,
        and sets the significance level threshold.

        Parameters:
        -----------
        bed_path : str
            File path to the .bed file.
        bim_path : str
            File path to the .bim file.
        fam_path : str
            File path to the .fam file.
        threshold : float
            Significance level threshold for statistical tests.
        """
        self.report = []
        self.bed_path = bed_path
        self.bim_path = bim_path
        self.fam_path = fam_path
        self.threshold = threshold
        print(f"association.py bed path -->{self.bed_path}")
        
        self.geno = IODataHandler._load_genotype_data(self)
        print(f"Load Genotype data-->{self.geno}")
        self.y = IODataHandler._load_phenotype_data(self)
        print(f"load Phenotype data-->{self.y}")
        self.snp_ids = self._load_snp_ids()
        print(f"association.py snp_ids -->{self.snp_ids}")
        self.X_normalized = self._preprocess_genotype_data()
        print(f"association.py normalized -->{self.X_normalized}")
        

    def log_report(self, message: str):
        '''helper function to generate report for the function executed'''
        print(message)
        self.report.append(message)
    def _load_snp_ids(self) -> np.ndarray:
        """
        Loads SNP IDs from the .bim file.

        Returns:
        --------
        np.ndarray:
            Array of SNP IDs.
        """
        self.log_report("Loading snp_ids")
        return pd.read_csv(self.bim_path, sep='\s+', header=None)[1].values

    def _preprocess_genotype_data(self) -> np.ndarray:
        """
        Handles missing values and normalizes genotype data.

        Returns:
        --------
        np.ndarray:
            Normalized genotype data as a NumPy array.
        """
        self.log_report("Preprocessing genotype data")
        imputer = SimpleImputer(strategy='mean')
        X_imputed = imputer.fit_transform(self.geno)

        scaler = StandardScaler()
        return scaler.fit_transform(X_imputed)

    def logistic_regression_comparison(self) -> pd.DataFrame:
        """
        Performs logistic regression on genotype data using both sklearn and statsmodels, 
        and compares the p-values for each SNP.

        Returns:
        --------
        pd.DataFrame:
            DataFrame containing SNP IDs and p-values from both sklearn and statsmodels.
        """
        p_values_sklearn = []
        p_values_statsmodels = []

        for snp_idx in range(self.X_normalized.shape[1]):
            snp_data = self.X_normalized[:, snp_idx].reshape(-1, 1)

            # Sklearn logistic regression
            p_value_sklearn = self._sklearn_logistic_regression(snp_data)
            p_values_sklearn.append(p_value_sklearn)

            # Statsmodels logistic regression
            p_value_statsmodels = self._statsmodels_logistic_regression(snp_data)
            p_values_statsmodels.append(p_value_statsmodels)
        self.log_report("Logistic regression Comparison function completed")
        return pd.DataFrame({
            'SNP ID': self.snp_ids,
            'p-value Statsmodels Logistic': p_values_statsmodels,
            'p-value Sklearn Logistic': p_values_sklearn
        })

    def linear_regression_comparison(self) -> pd.DataFrame:
        """
        Performs linear regression on genotype data using both sklearn and statsmodels,
        and compares the p-values for each SNP.

        Returns:
        --------
        pd.DataFrame:
            DataFrame containing SNP IDs and p-values from both sklearn and statsmodels.
        """
        p_values_sklearn = []
        p_values_statsmodels = []

        for snp_idx in range(self.X_normalized.shape[1]):
            snp_data = self.X_normalized[:, snp_idx].reshape(-1, 1)

            # Sklearn linear regression
            p_value_sklearn = self._sklearn_linear_regression(snp_data)
            p_values_sklearn.append(p_value_sklearn)

            # Statsmodels linear regression
            p_value_statsmodels = self._statsmodels_linear_regression(snp_data)
            p_values_statsmodels.append(p_value_statsmodels)
        self.log_report("Linear regression Comparison function completed")
        return pd.DataFrame({
            'SNP ID': self.snp_ids,
            'p-value Statsmodels Linear': p_values_statsmodels,
            'p-value Sklearn Linear': p_values_sklearn
        })

    def _sklearn_logistic_regression(self, snp_data: np.ndarray) -> float:
        """
        Performs logistic regression using sklearn and returns the p-value for the SNP.

        Parameters:
        -----------
        snp_data : np.ndarray
            Genotype data for a single SNP.

        Returns:
        --------
        float:
            p-value for the SNP coefficient.
        """
        try:
            model_sklearn = LogisticRegression(solver='liblinear')
            model_sklearn.fit(snp_data, self.y)

            coef = model_sklearn.coef_[0][0]
            pred_probs = model_sklearn.predict_proba(snp_data)[:, 1]

            X_design = np.hstack([np.ones((snp_data.shape[0], 1)), snp_data])
            V = np.diagflat(pred_probs * (1 - pred_probs))
            cov_matrix = np.linalg.inv(X_design.T @ V @ X_design)
            se = np.sqrt(np.diag(cov_matrix))
            z_scores = coef / se[1]
            #self.log_report("sklearn logistic regression Comparison function completed")
            return 2 * (1 - stats.norm.cdf(np.abs(z_scores)))
        except Exception:
            return np.nan

    def _statsmodels_logistic_regression(self, snp_data: np.ndarray) -> float:
        """
        Performs logistic regression using statsmodels and returns the p-value for the SNP.

        Parameters:
        -----------
        snp_data : np.ndarray
            Genotype data for a single SNP.

        Returns:
        --------
        float:
            p-value for the SNP coefficient.
        """
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=ConvergenceWarning)
                model_statsmodels = sm.Logit(self.y, sm.add_constant(snp_data)).fit(disp=False)
                #self.log_report("statsmodel logistic regression Comparison function completed")

            return model_statsmodels.pvalues[1]
        except Exception:
            return np.nan

    def _sklearn_linear_regression(self, snp_data: np.ndarray) -> float:
        """
        Performs linear regression using sklearn and returns the p-value for the SNP.

        Parameters:
        -----------
        snp_data : np.ndarray
            Genotype data for a single SNP.

        Returns:
        --------
        float:
            p-value for the SNP coefficient.
        """
        try:
            model_sklearn = LinearRegression()
            model_sklearn.fit(snp_data, self.y)
            coef = model_sklearn.coef_[0]

            y_pred = model_sklearn.predict(snp_data)
            residuals = self.y - y_pred
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((self.y - np.mean(self.y)) ** 2)
            r2 = 1 - (ss_res / ss_tot)
            if np.var(snp_data) == 0:
                return np.nan
            var_b = np.var(residuals) / (np.var(snp_data) * len(self.y))
            se_b = np.sqrt(var_b)

            if se_b == 0:
                return np.nan
            else:
                t_stat = coef / se_b
                return 2 * (1 - stats.t.cdf(np.abs(t_stat), df=len(self.y) - 2))
            #self.log_report("sklear logistic regression Comparison function completed")

        except Exception:
            return np.nan

    def _statsmodels_linear_regression(self, snp_data: np.ndarray) -> float:
        """
        Performs linear regression using statsmodels and returns the p-value for the SNP.

        Parameters:
        -----------
        snp_data : np.ndarray
            Genotype data for a single SNP.

        Returns:
        --------
        float:
            p-value for the SNP coefficient.
        """
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=ConvergenceWarning)
            model_statsmodels = sm.OLS(self.y, sm.add_constant(snp_data)).fit()
            #self.log_report("statsmodel linear regression Comparison function completed")

            return model_statsmodels.pvalues[1]
        except Exception:
            return np.nan

    def save_results(self, logistic_df: pd.DataFrame, linear_df: pd.DataFrame) -> None:
        """
        Saves the results of logistic and linear regression comparisons to a CSV file,
        including significance information based on Bonferroni correction.

        Parameters:
        -----------
        logistic_df : pd.DataFrame
            DataFrame containing logistic regression comparison results.
        linear_df : pd.DataFrame
            DataFrame containing linear regression comparison results.
        """
        comparison_df = pd.merge(logistic_df, linear_df, on='SNP ID')

        num_tests = len(self.snp_ids)
        alpha_bonferroni = self.threshold / num_tests

        comparison_df['Significant Statsmodels Logistic'] = comparison_df['p-value Statsmodels Logistic'] < alpha_bonferroni
        comparison_df['Significant Sklearn Logistic'] = comparison_df['p-value Sklearn Logistic'] < alpha_bonferroni
        comparison_df['Significant Statsmodels Linear'] = comparison_df['p-value Statsmodels Linear'] < alpha_bonferroni
        comparison_df['Significant Sklearn Linear'] = comparison_df['p-value Sklearn Linear'] < alpha_bonferroni

        comparison_df.to_csv('comparison_regression_results_with_significance.csv', index=False)
        self.log_report("Result is saved to a csv files")

        print(f"Comparison results with significance saved to 'comparison_regression_results_with_significance.csv'.")

    def generate_report(self, selected_functions):

        for func_name in selected_functions:
            func=getattr(self, func_name)
            if callable(func):
                if func_name =="logistic_regression_comparison":
                    df=func()
                    print(f"Logistic regression Comparison result -->{df.head()}")
                elif func_name=="linear_regression_comparison":
                    df_reg = func()
                    print(f"Linear Regression Comparison -->{df_reg.head()}")
                elif func_name=="_sklearn_logistic_regression":
                    p_values_sklearn = []
                    for snp_idx in range(self.X_normalized.shape[1]):
                        snp_data = self.X_normalized[:, snp_idx].reshape(-1, 1)

                        # Sklearn logistic regression
                        p_value_sklearn = func(snp_data)
                        p_values_sklearn.append(p_value_sklearn)
                    print(f"sklearn logistic regression --->", pd.DataFrame({
                    'SNP ID': self.snp_ids,
                    'p-value Sklearn Logistic': p_values_sklearn
                    }).head())
                    return pd.DataFrame({
                    'SNP ID': self.snp_ids,
                    'p-value Sklearn Logistic': p_values_sklearn
                    })
                elif func_name=="_statsmodels_logistic_regression":
                    p_values_statsmodels = []
                    for snp_idx in range(self.X_normalized.shape[1]):
                        snp_data = self.X_normalized[:, snp_idx].reshape(-1, 1)

                        p_value_statsmodels = func(snp_data)
                        p_values_statsmodels.append(p_value_statsmodels)
                    print(f"statsmodel_logistic_regression -->{pd.DataFrame({
                    'SNP ID': self.snp_ids,
                    'p-value Statsmodels Logistic': p_values_statsmodels,
                    }).head()}")
                    return pd.DataFrame({
                    'SNP ID': self.snp_ids,
                    'p-value Statsmodels Logistic': p_values_statsmodels,
                    })
                elif func_name =="_sklearn_linear_regression":
                    p_values_sklearn = []

                    for snp_idx in range(self.X_normalized.shape[1]):
                        snp_data = self.X_normalized[:, snp_idx].reshape(-1, 1)

                        # Sklearn linear regression
                        p_value_sklearn = self.func(snp_data)
                        p_values_sklearn.append(p_value_sklearn)
                    print(f"sklear_linear regression -->",pd.DataFrame({
                        'SNP ID': self.snp_ids,
                        'p-value Sklearn Linear': p_values_sklearn
                    }).head())
                    return pd.DataFrame({
                        'SNP ID': self.snp_ids,
                        'p-value Sklearn Linear': p_values_sklearn
                    })
                elif func_name =="_statsmodels_linear_regression":
                    p_values_statsmodels = []

                    for snp_idx in range(self.X_normalized.shape[1]):
                        snp_data = self.X_normalized[:, snp_idx].reshape(-1, 1)

                        # Statsmodels linear regression
                        p_value_statsmodels = self.func(snp_data)
                        p_values_statsmodels.append(p_value_statsmodels)
                    print(f"statsmodel_linear_regression --->",pd.DataFrame({
                        'SNP ID': self.snp_ids,
                        'p-value Statsmodels Linear': p_values_statsmodels
                    }).head())
                    return pd.DataFrame({
                        'SNP ID': self.snp_ids,
                        'p-value Statsmodels Linear': p_values_statsmodels
                    })
def main():
    avialable_functions=[
        'logistic_regression_comparison',
        'linear_regression_comparison',
        '_sklearn_logistic_regression',
        '_statsmodels_logistic_regression',
        '_sklearn_linear_regression',
        '_statsmodels_linear_regression'

    ]

    for idx, func in enumerate(avialable_functions,1):
        print(f"{idx},  {func}")
    selected_indices = input("Enter the number of the functions to run e.g 1 3 for the first and third function ").split()

    try:
        selected_functions=[avialable_functions[int(i)-1] for i  in selected_indices]
    except(ValueError, IndexError):
        print("Invalid selection, Please enter valid numbers")
        return
    
    bed_path="C:/Users/smith/OneDrive/Documents/GitHub/Fed-GWAS/data/begin.cc"
    bim_path="C:/Users/smith/OneDrive/Documents/GitHub/Fed-GWAS/data/begin.cc.bim"
    fam_path="C:/Users/smith/OneDrive/Documents/GitHub/Fed-GWAS/data/begin.cc.fam"
    threshold=0.05

    association = GWASAssociation(bed_path=bed_path, bim_path=bim_path,fam_path=fam_path, threshold=threshold)
    association.generate_report(selected_functions)
if __name__ == "__main__":
    main()
