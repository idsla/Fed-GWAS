from fedgwas.io import reader
from fedgwas.quality_control import qc
from pysnptools.snpreader import Bed
import pandas as pd


def main():
    print("hello World")


if __name__== "__main__":
    main()
    filename = 'data/begin.cc'
    bed=reader.snp_reader_Bed(filename)
    num_individuals=bed.iid_count
    num_snps=bed.sid_count
    #print(bed.get_fam())
    #qc= qc.QualityControl()
    #df=qc.calculate_maf(bed)
    #bim, fam, G = reader.pandas_plink(filename)
    #bim_df=pd.DataFrame(bed.get_bim()).reset_index()
    fraction = 0.2
    print(f"Number of individuals :{num_individuals} and num of snps: {num_snps}")
    num_snps = int(num_snps*0.2)
    num_individuals = int(num_individuals*0.2)
    print(f"Number of individuals :{num_individuals} and num of snps: {num_snps}")
    subset = bed[:num_snps, num_individuals].read()
    print(f"G subset {subset}")
    data_A = subset[:, :num_individuals // 2].val.T
    data_B = subset[:, num_individuals // 2:].val.T
    print(f"Data A {data_A} and Data B {data_B}")


    