QUALITY_CONTROL = {
    'filter_missingness_samples':{
        'threshold': 0.05,     # threshold for filtering missingness
    },
    'hardy_weinberg_test':{
        'threshold': 1e-6 ,     # threshold for HWE 
    },
    'geno':{
        'threshold': 0.05,     # threshold for filtering missingness
    }
}

KINSHIP = {
    "final_kinship_estimate": {
        "firstkin": 0.35,   # threshold for first degree kinship
        'secondkin': 0.18,  # threshold for second degree kinship
        'thirdkin': 0.09,   # threshold for third degree kinship
    }
}

ASSOCIATION = {

}