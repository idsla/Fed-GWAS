from pyplink import PyPlink
from pandas_plink import read_plink
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed

filename = 'data/begin.cc.bed'

def read_file(filename):
    '''
    Returns bed_file
    '''
    bed_file = PyPlink(filename)
    return bed_file
def pandas_plink(filename):
    """
    Reads PLINK genotype data from a specified file and returns the BIM, FAM, and genotype matrices.

    Parameters:
    filename (str): The path to the PLINK file prefix (without the .bed, .bim, or .fam extensions).

    Returns:
    tuple: A tuple containing three pandas DataFrames - (bim, fam, G).
        bim (DataFrame): SNP information from the BIM file.
        fam (DataFrame): Individual information from the FAM file.
        G (DataFrame): Genotype matrix.
    """
    (bim, fam, G) = read_plink(filename)
    return bim, fam, G
def bim_to_bim_df(bedfile):
    '''
    Returns bim Data Frame from .bed file
    '''
    bim_df = pd.DataFrame(bedfile.get_bim()).reset_index()
    return bim_df

def fam_to_fam_df(bedfile):
    '''
    Returns fam Data Frame from .bed file
    '''
    fam_df = pd.DataFrame(bedfile.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
    return fam_df

def snp_reader_Bed(filepath):
    snp_reader= Bed(filepath)
    return snp_reader 
