from fedgwas.association.association import GWASAssociation
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
# if __name__ == "__main__":
#     assoc = GWASAssociation(
#         bed_path='/data/begin.cc',
#         threshold=0.05

#     )
#     logistic_df = assoc.logistic_regression_comparison()
#     linear_df = assoc.linear_regression_comparison()
#     assoc.save_results(logistic_df, linear_df)