from fedgwas.association.association import GWASAssociation

if __name__ == "__main__":
    assoc = GWASAssociation(
        bed_path='/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc6.bed',
        bim_path='/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc6.bim',
        fam_path='/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc6.fam',
        threshold=0.05

    )
    logistic_df = assoc.logistic_regression_comparison()
    linear_df = assoc.linear_regression_comparison()
    assoc.save_results(logistic_df, linear_df)