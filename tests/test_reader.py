import logging as logging
import datetime
import pandas as pd
from  fedgwas.quality_control.qc import QualityControl 
from fedgwas.io.reader import IODataHandler
from fedgwas.quality_control.qc_utils import QCUtils
from fedgwas.parameters import QUALITY_CONTROL

#Setup logging
logging.basicConfig(filename='reader_IO.log', level=logging.INFO, format='%(message)s')
start_time = datetime.datetime.now()
logging.info(f"Start time: {start_time.strftime('%a %b %d %H:%M:%S %Y')}")

if __name__ == "__main__":
    # Setup argument parser

        io = IODataHandler(
            bed_path='/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc6.bed',
            bim_path='/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc6.bim',
            fam_path='/Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc6.fam'
        )

        # Read the .fam and .bim files
        fam_df = io.read_fam()
        bim_df = io.read_bim()

        # Read the .bed file
        genotype_data, bed = io.read_bed()
        geno_df = pd.DataFrame(genotype_data)
        # Print the results
        logging.info(f"Genotype data:\n{genotype_data}")
        qcInstance = QualityControl()
        missing_threshold = QUALITY_CONTROL['qc_threshold']['missing_threshold']
        hwe_threshold = QUALITY_CONTROL['qc_threshold']['hwe_threshold']
        qcInstance.calculate_missing_rate(geno_df, bed, 'missing_rate')
        qcInstance.filter_missingness_samples(geno_df, fam_df, missing_threshold, 'individual_missing')
        qcInstance.geno(genotype_data, bim_df, missing_threshold, 'snp_missing')
        filtered_geno, snp_indices_to_keep = qcInstance.hardy_weinberg_test(geno_df, hwe_threshold)
        io.save_filtered_data(bim_df, filtered_geno, snp_indices_to_keep, 'hwe_test')