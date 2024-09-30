
from fedgwas.quality_control.qc import QualityControl
def main():

    available_functions=[
        'check_sex',
        'filter_maf_variants',
        'calculate_missing_rate',
        'hardy_weinberg_test',
        'filter_missingness_samples',
        'geno'

    ]

    for idx, func in enumerate(available_functions,1):
        print(f"{idx}.  {func}")

    selected_indices= input("Enter the numbers of the functions to run (e.g. 1 3 for this first and third) ").split()

    try: 
        #covert selected indices to function names
        selected_functions=[available_functions[int(i)-1] for i in selected_indices]
    except(ValueError, IndexError):
        print("Invalid Selection, please enter valid numbers")
        return

    bed_path="C:/Users/smith/OneDrive/Documents/GitHub/Fed-GWAS/data/begin.cc"
    bim_path="data/begin.cc.bim"
    fam_path="data/begin.cc.fam"
    threshold=0.05
    
    qc = QualityControl(bed_path=bed_path, bim_path=bim_path, fam_path=fam_path, threshold=threshold)

    qc.run_quality_control_pipeline(selected_functions)

if __name__=="__main__":
    main()
# import pytest
# import pandas as pd
# from fedgwas.quality_control.qc import QualityControl

# @pytest.fixture
# def genotype_data():
#     data = None  # Replace this with your sample data
#     return pd.DataFrame(data)

# @pytest.fixture
# def phenotype_data():
#     data = None  # Replace this with your sample data
#     return pd.DataFrame(data)

# def test_check_sex(genotype_data, phenotype_data):
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_calculate_missing_rate(genotype_data):
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_calculate_maf(genotype_data):
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_hardy_weinberg_test(genotype_data):
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_filter_missingness_samples(genotype_data):
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_filter_maf_variants(genotype_data):
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_filter_data():
#     qc = QualityControl()
#     pass # Replace this with your test code

# def test_generate_report():
#     qc = QualityControl()
#     pass # Replace this with your test code
