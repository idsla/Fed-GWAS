from fedgwas.io import reader
from fedgwas.quality_control import qc
from pysnptools.snpreader import Bed
import pandas as pd
import numpy as np
from pyplink import PyPlink


def main():
    print("hello World")


if __name__== "__main__":
    main()
    filename = 'data/begin.cc'
    bed=reader.snp_reader_Bed(filename)
    bedpy= PyPlink(filename)
    bimdf= bedpy.get_bim()
    print(bimdf)
    chrom=bimdf['chrom']
    a1=bimdf['a1']
    a2=bimdf['a2']
    pos=bimdf['pos']
    snps, mafs, nchrobs = [], [], []
    
    for i in range(bed.sid_count):
        # Read the genotype data for the SNP
        genotype = bed[:, i].read().val.flatten()
        
        # Calculate allele counts
        count_a1 = 2 * (genotype == 0).sum() + (genotype == 1).sum()
        count_a2 = 2 * (genotype == 2).sum() + (genotype == 1).sum()
        no_of_chzms = count_a1 + count_a2
        
        # Calculate allele frequencies
        freq_a1 = count_a1 / no_of_chzms
        freq_a2 = count_a2 / no_of_chzms
        
        # Calculate MAF
        maf_coeff = min(freq_a1, freq_a2)
        mafs.append(maf_coeff)
        nchrobs.append(no_of_chzms)
        snps.append(bed.sid[i])

    df = pd.DataFrame({
        'CHR': chrom,
        'SNP': snps,
        'POS': pos,
        'A1': a1,
        'A2': a2,
        'MAF': mafs,
        'NCHROBS': nchrobs
    })
    print(df.head())
    # bedpy= PyPlink(filename)
    # print(bedpy.get_bim())

    # #print(f"Genotype array---> {bed.read().val}")
    # genotype_arr= bed.read().val
    # nan_count= np.sum(np.isnan(genotype_arr))
    # genotype_array = np.nan_to_num(genotype_arr, nan=-1).astype(int)
    # #print(f"nan count {nan_count}")
    # #print(np.unique(genotype_array))
    # genotype_data = bed.read()
    # bed_file=reader.snp_reader_Bed(filename)
    # # Check the length of each column
    # print(len(bed_file.pos[:, 0]))  # Chromosome
    # print(len(bed_file.sid))        # SNP ID
    # print(len(bed_file.pos[:, 1]))  # Genetic distance
    # print(len(bed_file.pos[:, 2]))  # Position in base pairs
    # print(len(bed_file.row[:, 0]))  # Allele 1
    # print(len(bed_file.row[:, 1]))  # Allele 2

    # print(f"Chromosome bed_file.pos[:, 0] ---> {bed_file.pos[:, 0]} len --> {len(bed_file.pos[:, 0])}")  # Chromosome
    # print(f"SNP ID bed_file.sid --->{bed_file.sid} len--->{len(bed_file.sid)}")        # SNP ID
    # print(f"Genetic distance bed_file.pos[:, 1] -->{bed_file.pos[:, 1]} len---> {len(bed_file.pos[:, 1])}")  # Genetic distance
    # print(f"Position in base pairs  (bed_file.pos[:, 2])--->{(bed_file.pos[:, 2])} len---> {len(bed_file.pos[:, 2])}")  # Position in base pairs
    # print(f"Allele 1 bed_file.row[:, 0]---->{bed_file.row[:, 0]} len--->{len(bed_file.row[:, 0])}")  # Allele 1
    # print(f"Allele 2 bed_file.row[:, 1] ---->{bed_file.row[:, 1]} len--->{len(bed_file.row[:, 1])}")  # Allele 2

    # Check the data to ensure it's loaded correctly
    #print(bim_df.head())
    # print(f"chrom--->{genotype_data.pos[:,0]} --->length : {len(genotype_data.pos[:,0])}")
    # print(f"snp --->{genotype_data.sid} ---> length : {len(genotype_data.sid)}")
    # print(f"genetic dist ---> {genotype_data.pos[:,1]}  ")
    # print(f"pos ----> {genotype_data.pos[:,2]} length----->:{len(genotype_data.pos[:,2])}")
    # print(f"a1----> {genotype_data.row[:,0]}, length ----> {len(genotype_data.row[:,0])}")
    # print(f"a2 ----> {genotype_data.row[:,1]}, length ----> {len(genotype_data.row[:,1])}")
    # bim_df= pd.DataFrame({
    #     'chrom': genotype_data.pos[:,0],
    #     'snp': genotype_data.sid,
    #     'genetic dist': genotype_data[:,1],
    #     'pos':genotype_data.pos[:,2],
    #     'a1':genotype_data[:,0],
    #     'a2':genotype_data.row[:,1]
    # })
    # maf_values=[]
    # snps=[]
    # nchrobs=[]
    # #allele_counts= np.apply_along_axis(lambda x: np.bincount(x, minlength=3), axis=0, arr= genotype_array)
    # for snp in range(genotype_array.shape[1]):
    #     valid_genotypes = genotype_array[:, snp][genotype_array[:, snp]>=0]

    #     if valid_genotypes.size >0:
    #         no_of_chzms = 2 * np.sum(allele_counts)
    #         allele_counts= np.bincount(valid_genotypes, minlength=3)
    #         n_allele_counts= 2*np.sum(allele_counts)
    #         p_freq = (2*allele_counts[2] + allele_counts[1])/n_allele_counts
    #         q_freq = 1-p_freq

    #         maf = np.minimum(p_freq, q_freq)
    #     else:
    #         maf = np.nan
    #     maf_values.append(maf)
    #     nchrobs.append(no_of_chzms)
    #     snps.append(bim_df['snp'][snp])

    # print(f" Length of snps is{len(snps)} --> length of nchrobs is {len(nchrobs)} --> length of snps is {len(snps)}")
#     print(f"maf values--->{maf_values}")


#     df = pd.DataFrame({
#     'CHR': bim_df['chrom'],
#     'SNP': snps,
#     'POS': bim_df['pos'],
#     'A1': bim_df['a1'],
#     'A2': bim_df['a2'],
#     'MAF': maf_values,
#     'NCHROBS': nchrobs
# })
# Save the DataFrame to a CSV file
#df.to_csv('snp_maf_data.csv', index=False)   
   # print(f"Allele counts-->{allele_counts}")
    # num_individuals=bed.iid_count
    # num_snps=bed.sid_count
    #print(bed.get_fam())
    #qc= qc.QualityControl()
    #df=qc.calculate_maf(bed)
    #bim, fam, G = reader.pandas_plink(filename)
    #bim_df=pd.DataFrame(bed.get_bim()).reset_index()
    # fraction = 0.2
    # print(f"Number of individuals :{num_individuals} and num of snps: {num_snps}")
    # num_snps = int(num_snps*0.2)
    # num_individuals = int(num_individuals*0.2)
    # print(f"Number of individuals :{num_individuals} and num of snps: {num_snps}")
    # subset = bed[:num_snps, num_individuals].read()
    # print(f"G subset {subset}")
    # data_A = subset[:, :num_individuals // 2].val.T
    # data_B = subset[:, num_individuals // 2:].val.T
    # print(f"Data A {data_A} and Data B {data_B}")


    