import pandas as pd
import numpy as np
import os

# filter functions
def snp_genotype_filter(bfile: str, threshold: float, outfile: str):
    command = "./plink --bfile {} --geno {} --silent --make-bed --out {}".format(bfile, threshold, outfile)
    os.system(command)

def ind_genotype_filter(bfile: str, cutoff: float, outfile: str):
    command = "./plink --bfile {} --mind {} --silent --make-bed --out {}".format(bfile, cutoff, outfile)
    os.system(command)

def individual_filter(bfile: str, keepfile: str, outfile: str):
    #--keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis
    command = "./plink --bfile {} --keep {} --silent --make-bed --out {}".format(bfile, keepfile, outfile)
    os.system(command)

def impute_sex(bfile: str):
    outfile = bfile + "_sex_imputed"
    command = "./plink --bfile {} --impute-sex --silent --make-bed --out {}".format(bfile, outfile)
    os.system(command)

def remove_sex(bfile: str, remove_file: str, out: str):
    if os.path.isfile(remove_file):
        command = "./plink --bfile {} --remove {} --make-bed --out {}".format(bfile, remove_file, out)
        os.system(command)
    else:
        print("error in remove_sex function! {} file is not found!").format(remove_file)

def autosomal_snps_filter(bfile: str, outfile: str):
    bim_file = bfile + ".bim"
    bim = pd.read_csv(bim_file, delimiter="\t", skipinitialspace=True)
    bim.columns = ['chrom', 'snp', 'cm', 'pos', 'a0', 'a1']
    bim['chrom'] = bim['chrom'].astype('int32')
    bim_filtered = bim.loc[(bim['chrom'] >= 1) & (bim['chrom'] <= 22)]
    bim_filtered['snp'].to_csv(outfile, index=None, sep=' ')
    return bim_filtered

def select_autosomal_snp(bfile: str, auto_file: str, outfile: str):
    # takes in output from "automsomal_snps_filter"
    command = "./plink --bfile {} --extract {} --silent --make-bed --out {}".format(bfile, auto_file, outfile)
    os.system(command)

def maf_filter(bfile: str, threshold: float, outfile: str):
    command = "./plink --bfile {} --maf {} --make-bed --out {}".format(bfile, threshold, outfile)
    os.system(command)

def hardy_weinberg_filter(bfile: str, outfile: str, control: bool, threshold: float):
    if control:
        command = "./plink --bfile {} --hwe {} --silent --make-bed --out {}".format(bfile, threshold, outfile)
    else:
        command = "./plink --bfile {} --hwe include-nonctrl {} --silent --make-bed --out {}".format(bfile, threshold, outfile)
    os.system(command)
    #os.system(command2)

def ld_pruning_filter(bfile: str, snp_file: str, outfile: str, window: int=50, shift: int=5, correlation_threshold: int=0.2, correlation_method: str="pairwise"):
    if correlation_method=="multiple":
        # add check for correlation threshold if "multiple" is chosen - must be > 1
        command = "./plink --bfile {} --indep {} {} {} --out {}".format(bfile, window, shift, correlation_threshold, snp_file)
    elif correlation_method=="pairwise":
        command = "./plink --bfile {} --indep-pairwise {} {} {} --out {}".format(bfile, window, shift, correlation_threshold, snp_file)
    snp_in = snp_file + ".prune.in"
    command2 = "./plink --bfile {} --extract {} --het --out {}".format(bfile, snp_in, outfile)
    os.system(command)
    os.system(command2)

def heterozygosity_filter(bfile: str, failed_file: str, outfile: str):
    command = "./plink --bfile {} --remove {} --make-bed --out {}".format(bfile, failed_file, outfile)
    os.system(command)

def relatedness_filter(bfile: str, remove_file: str, outfile: str):
    command = "./plink --bfile {} --remove {} --make-bed --out {}".format(bfile, remove_file, outfile)
    os.system(command)

def founders_filter(bfile: str, outfile: str):
    command = "./plink --bfile {} --filter-founders --silent --make-bed --out {}".format(bfile, outfile)
    os.system(command)
    return outfile

def sample_ids_filter(df, column, filtered):
    '''Return sample IDs from one dataframe (df) that are present in another dataframe (filtered).
    :param df: dataframe that contains the samples to be filtered
    :param column: column that contains the sample IDs
    :param filtered: dataframe that contains the filtered samples
    :return: Pandas DataFrame object
    :rtype object
    >>> get_sample_ids(original_df, 'iid', fitlered_df)
    ids'''
    ids = df[~df[column].isin(filtered[column].tolist())][column]
    return ids

def rename_filter(bfile: str, outfile: str):
    command = "./plink --bfile {} --silent --make-bed --out {}".format(bfile, outfile)
    os.command(command)
