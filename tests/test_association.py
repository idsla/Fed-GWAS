from fedgwas.association.association import GWASAssociation

if __name__ == "__main__":
    assoc = GWASAssociation(
        bed_path='/data/begin.cc',
        threshold=0.05

    )
    logistic_df = assoc.logistic_regression_comparison()
    linear_df = assoc.linear_regression_comparison()
    assoc.save_results(logistic_df, linear_df)