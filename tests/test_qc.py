import argparse
import pandas as pd
import numpy as np
import logging
from pyplink import PyPlink
import sys
from fedgwas.quality_control.qc import QualityControl
#from fedgwas.parameters import QUALITY_CONTROL
sys.path.insert(0,'C://Users//smith//OneDrive//Documents//GitHub//Fed-GWAS//fedgwas')
logging.basicConfig(filename="qctest.log", filemode='w', level= logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def main():

    parser= argparse.ArgumentParser( description="Run QC from command line for Test")

    parser.add_argument(
        '--function',
        type=str,
        required=True,
        choices=['check_sex', 'calculate_missing_rate', 'calculate_maf'],
        help="command line commands"
    )

    parser.add_argument('--bed', type= str, help="Path of the file")
    args= parser.parse_args()
    qc= QualityControl()
    bed= None
    if args.bed:
        bed= PyPlink(args.bed)

    if args.function =='check_sex':
        print("Running check SEX function")
        qc.check_sex(bed)
    elif args.function == 'calculate_maf':
        print("Running Calculate MAF Rate function")
        qc.calculate_maf(bed)

if __name__=="__main__":
    main()