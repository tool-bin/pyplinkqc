import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib.backends.backend_pdf import PdfPages

# analysis functions
def calculate_missingness(df, column, threshold):
    '''Identify SNPs with a missingness rate below a certain threshold.
    :param df: dataframe containing all SNPs and their missingness rate
    :param column: column that contains missingness rate values within df
    :param threshold: threshold for removing SNPs above missingness rate
    :return: Pandas DataFrame object
    :rtype object
    >>> calculate_missingess(missigness_df, 'F_MISS', 0.05)
    missing'''
    missing = df.loc[df[column] < threshold]
    print("total removed: ", df.shape[0] - missing.shape[0])
    return missing

def get_sample_ids(df, column, filtered):
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

def heterozygosity_ind_report(infile: str, outfile: str):
    #infile should be the het_check.het file from the heterozygosity_report function
    het_check = pd.read_csv(infile, delimiter=" ", skipinitialspace=True)
    het_check['het_rate'] = (het_check['N(NM)'] - het_check['O(HOM)']) / het_check['N(NM)']
    het_check['low_limit'] = het_check['het_rate'].mean() - (3 * het_check['het_rate'].std())
    het_check['up_limit'] = het_check['het_rate'].mean() + (3 * het_check['het_rate'].std())
    het_check['outlier'] = np.where((het_check['het_rate'] > het_check['up_limit']) |  (het_check['het_rate'] < het_check['low_limit']),1, 0)
    het_fail = het_check.loc[het_check['outlier'] == 1]
    #print("original number of samples: {}".format(het_check.shape[0]))
    print("number of samples with deviating heterozygosity rates: {}".format(het_fail.shape[0]))

    het_fail.iloc[:, :2].to_csv(outfile, index=None, sep=' ')
    return het_check

def missingness_report(bfile: str, outfile: str):
    command = "./plink --bfile {} --missing --silent --out {}".format(bfile, outfile)
    os.system(command)

def check_sex(bfile: str):
    command = "./plink --bfile {} --check-sex --silent".format(bfile)
    os.system(command)

def check_sex_report(file: str="plink.sexcheck"):
    sexcheck = pd.read_csv(file, delimiter=" ", skipinitialspace=True)
    problems = sexcheck.loc[sexcheck['STATUS'] == "PROBLEM"]
    print("number of sex discrepancies: ", problems.shape[0])
    problems[['FID', 'IID']].to_csv("sex_discrepancy.txt", index=None, sep=' ')
    return problems

def relatedness_report(bfile, indep_snp_file, threshold=0.2):
    command = "./plink --bfile {} --extract {} --genome --min {} --silent --memory --out pihat_min_{}".format(bfile, indep_snp_file, threshold, threshold)
    os.system(command)

def maf_check_report(bfile):
    command = "./plink --bfile {} --freq --silent --out MAF_check".format(bfile)
    os.system(command)

def heterozygosity_report(bfile, outfile):
    #bfile should be ld_check based on output of ld_pruning_filter function
    #this is because the heterozygosity rates should be calculated from uncorrelated regions of the genome, so we exclude retions of high LD
    command = "./plink --bfile {}  --het --silent --out {}".format(bfile, outfile)
    os.system(command)

def hardy_weinberg_report(bfile):
    command = "./plink --bfile {} --hardy".format(bfile)
    os.system(command)

def relatedness_check_report(bfile: str, indep_snp_file: str, outfile: str, threshold: float=0.2):
    command = "./plink --bfile {} --extract {} --silent --genome --min {} --out {}".format(bfile, indep_snp_file, threshold, outfile)
    os.system(command)

def relatives_low_call_rate_report(imissfile: str, relatfile: str, outfile: str):
    imissfile = imissfile + ".imiss"
    imiss = pd.read_csv(imissfile, delimiter=" ", skipinitialspace=True)
    relat = pd.read_csv(relatfile, delimiter=" ", skipinitialspace=True)
    fid_first_list = relat['IID1'].tolist()
    fid_second_list = relat['IID2'].tolist()
    relat_miss1 = relat.set_index('IID1').join(imiss.set_index('IID')[['F_MISS']])
    relat_miss2 = relat_miss1.rename_axis('IID1').reset_index().set_index('IID2').join(imiss.set_index('IID')[['F_MISS']], rsuffix="_2").rename_axis('IID2').reset_index()
    ids = np.where(relat_miss2['F_MISS']>relat_miss2['F_MISS_2'], relat_miss2['IID1'], relat_miss2['IID2'])
    # first_selected = imiss[imiss['IID'].isin(fid_first_list)].reset_index(drop=True)
    # second_selected = imiss[imiss['IID'].isin(fid_second_list)].reset_index(drop=True)
    # selected = np.where(first_selected['F_MISS'] > second_selected['F_MISS'],first_selected[['FID', 'IID']],
    #                 second_selected[['FID', 'IID']])
    selected = imiss.loc[imiss['IID'].isin(ids)][['FID', 'IID']]
    np.savetxt(outfile, selected, delimiter=' ',  fmt='%d %s')
    return selected

def write_fail_file(ids_failed, filename="failed_ids"):
    '''Write either failed sample IDs or SNPs to two files.
    One is a human-readable csv file, the other is a csv file intended to be read into a Pandas DataFrame.
    :param ids_failed: dictionary that contains the QC checks and the corresponding failed ids.
    :param filename: name of the file to be written to
    >>> write_fail_file(ids_failed, failed_ids)
    '''
    hr = filename + "_hr.csv"
    pr = filename + "_pr.csv"
    with open(hr, "w") as f:
        for k,v in ids_failed.items():
            f.write(str(k) + ": " + str(v) + "\n")
    pd.DataFrame.from_dict(data=ids_failed, orient='index').to_csv(pr, header=True)

def save_pdf(infile, plots):
    ''' Write a list of matplotlib Figure objects to a pdf file in the current working directory.
    :param infile: prefix of the plink input files
    :param plots: list of Figure objects to written to pdf file
    >>> save_pdf("HapMap_3_r3_1", figures)
    "HapMap_3_r3_1_qc_report.pdf"'''

    title = infile + "_report.pdf"

    with PdfPages(title) as pdf:
        for plot in plots:
            pdf.savefig(plot)

def snps_failed_report(write: bool=False, miss_threshold: float=0.2, maf_threshold: float=0.00001, hwe_threshold: float=1e-6, lmiss_file: str="plink.lmiss", maf_file: str="MAF_check.frq", hwe_file: str="plink.hwe"):
    snps = {}
    ids_list = []
    lmiss = pd.read_csv(lmiss_file, delimiter=" ", skipinitialspace=True)

    missing_snps = lmiss.loc[lmiss['F_MISS'] > miss_threshold]
    snps['missing_snps'] = missing_snps['SNP'].tolist()
    ids_list.append(missing_snps['SNP'].tolist())
    print("total missing snps failed: ", len(missing_snps['SNP'].tolist()))

    # MAF
    maf = pd.read_csv(maf_file, delimiter=" ", skipinitialspace=True)
    rare = maf.loc[maf['MAF'] < maf_threshold]
    snps['maf'] = rare['SNP'].tolist()
    ids_list.append(rare['SNP'].tolist())
    print("total maf snps failed: ", len(rare['SNP'].tolist()))

    # HWE departures
    hardy = pd.read_csv(hwe_file, delimiter=" ", skipinitialspace=True)
    hwe_failed = hardy.loc[hardy['P'] < hwe_threshold]
    snps['hwe'] = hwe_failed['SNP'].tolist()
    ids_list.append(hwe_failed['SNP'].tolist())
    print("total hwe snps failed: ", len(hwe_failed['SNP'].tolist()))

    # graph everything
    tests = ['SNP Missingness', 'Minor Allele Frequency', 'Outlying HWE']
    fail_counts = [len(missing_snps['SNP'].tolist()), len(rare['SNP'].tolist()), len(hwe_failed['SNP'].tolist())]
    total_fails = set(x for l in ids_list for x in l)
    print("total fails: ", len(total_fails))

    fig = plt.figure(figsize=(8,6))
    plt.tight_layout()
    plt.bar(x=tests, height=fail_counts)
    plt.title("SNPs failing QC checks (total: {}/{})".format(len(total_fails), lmiss.shape[0]))
    plt.xlabel("QC Test")
    plt.ylabel("Number of SNPs")
    plt.tick_params(axis='x', rotation=90)

    if write:
        write_fail_file(snps, "failed_snps_ids")

    return fig

def sample_failed_report(write=True, miss_threshold=0.2, imiss_file="plink.imiss", lmiss_file="plink.lmiss", sexcheck_file="plink.sexcheck", het_failed_file="heterozygosity_failed.txt", ibd_state=True, ibd_file="pihat_min0.2.genome"):
    ids = {}
    ids_list = []

    # SNP missingness
    imiss = pd.read_csv(imiss_file,  delimiter=" ", skipinitialspace=True)
    lmiss = pd.read_csv(lmiss_file,  delimiter=" ", skipinitialspace=True)
    ind_missing_filtered = calculate_missingness(imiss, 'F_MISS', miss_threshold)
    missing_ids = get_sample_ids(imiss, 'IID', ind_missing_filtered)

    ids['missing'] = missing_ids.tolist()
    ids_list.append(missing_ids.tolist())

    # mismatched sex
    sex = pd.read_csv(sexcheck_file, delimiter=" ", skipinitialspace=True)
    sex_mismatches = sex.loc[sex['STATUS'] == "PROBLEM"]
    sex_mismatches_counts = sex['STATUS'].value_counts()
    print("total sex mismatches: ", sex_mismatches.shape[0])
    sex_mismatches_ids = sex_mismatches['IID'].tolist()

    ids['sex_mismatches'] = sex_mismatches_ids
    ids_list.append(sex_mismatches_ids)

    # outlying heterozygosity
    het_failed = pd.read_csv(het_failed_file, delimiter=" ")
    het_failed_ids = het_failed['IID'].tolist()
    print("total het failed mismatches: ", len(het_failed_ids))

    ids['het_failed'] = het_failed_ids
    ids_list.append(het_failed_ids)

    if ibd_state:
        # high IBD - pi_hat threshold
        ibd = pd.read_csv(ibd_file, delimiter=" ", skipinitialspace=True)
        print("total ibd failures: ", ibd.shape[0])
        ibd_ids = ibd['IID1'].tolist()
        ids['relatedness_failed'] = ibd_ids
        ids_list.append(ibd_ids)
        # graph everything
        tests = ['SNP Missingness', 'Sex Mismatches', 'Outlying Heterozygosity', 'Cryptic Relatedness']
        fail_counts = [len(missing_ids), len(sex_mismatches_ids), len(het_failed_ids), len(ibd_ids)]
        total_fails = set(x for l in ids_list for x in l)
        print("total samples failed: {}/{}".format(len(total_fails), imiss.shape[0]))

    else:
        # graph everything
        tests = ['SNP Missingness', 'Sex Mismatches', 'Outlying Heterozygosity']
        fail_counts = [len(missing_ids), len(sex_mismatches_ids), len(het_failed_ids)]
        total_fails = set(x for l in ids_list for x in l)
        print("total samples failed: {}/{}".format(len(total_fails), imiss.shape[0]))

    fig = plt.figure(figsize=(8,6))
    plt.tight_layout()
    plt.bar(x=tests, height=fail_counts)
    plt.title("Samples failing QC checks (total: {}/{})".format(len(total_fails), imiss.shape[0]))
    plt.xlabel("QC Test")
    plt.ylabel("Number of samples")
    plt.tick_params(axis='x', rotation=90)

    if write:
        write_fail_file(ids, "failed_sample_ids")

    return fig
