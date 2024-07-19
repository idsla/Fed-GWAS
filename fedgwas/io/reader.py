from pyplink import PyPlink
import numpy as np
import pandas as pd

filename = 'C:/Users/smith/Downloads/handsonData/begin.cc'

def read_file(filename):
    bed_file = PyPlink(filename)
    return bed_file

def bim_to_bim_df(bedfile):
    bim_df = pd.DataFrame(bedfile.get_bim()).reset_index()
    return bim_df

def fam_to_fam_df(bedfile):
    fam_df = pd.DataFrame(bedfile.get_fam(), columns=['fid', 'iid', 'father', 'mother', 'gender', 'status'])
    return fam_df
