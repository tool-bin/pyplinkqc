import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib.backends.backend_pdf import PdfPages
import os

def plot_missingness_hist(missfile: str="plink"):
    """Plot histograms of SNP missingness for samples and SNPs.

    Input files should be generated by report.missingness() function.

    Key arguments:
    --------------
    missfile: str
        prefix for the plink file containing missingess information

    Returns:
    --------
    Figure object
    """
    imiss_file = missfile+".imiss"
    lmiss_file = missfile+".lmiss"
    imiss = pd.read_csv(imiss_file, delimiter=" ", skipinitialspace=True)
    lmiss = pd.read_csv(lmiss_file, delimiter=" ", skipinitialspace=True)

    fig, ax = plt.subplots(1,2, figsize=(16,6), sharex=True)
    imiss.hist(column='F_MISS', ax=ax[0])
    ax[0].axvline(0.2, c='red', linestyle='--')
    ax[0].set_xlabel("Proportion of missing SNPs")
    ax[0].set_ylabel("Number of individuals")
    ax[0].set_title("Proportion of missing SNPs per individual \n (> 0.2 are removed)")
    lmiss['F_MISS'].hist(ax=ax[1])
    ax[1].axvline(0.2, c='red', linestyle='--')
    ax[1].set_xlabel("Proportion of individuals with missing SNPs")
    ax[1].set_ylabel("Number of SNPs")
    ax[1].set_title("Proportion of missing individuals per SNP \n (> 0.2 are removed)")
    return fig

def plot_sexcheck_hist(sexcheckfile: str="plink.sexcheck"):
    """Plot histograms of inbreeding coefficents for reported females/males.

    Input files should be generated by qc_report.check_sex().

    Key arguments:
    --------------
    sexcheckfile: str
        prefix for the plink file containing sex check information

    Returns:
    --------
    Figure object
    """
    sexcheck = pd.read_csv(file, delimiter=" ", skipinitialspace=True)
    females = sexcheck.loc[sexcheck['PEDSEX'] == 2]
    males = sexcheck.loc[sexcheck['PEDSEX'] == 1]

    fig, ax = plt.subplots(1, 2, figsize=(16,6))
    females.hist(column='F', ax=ax[0])
    ax[0].axvline(0.2, c='r', ls='--')
    ax[0].set_ylabel("Number of Individuals")
    ax[0].set_xlabel("F Value (X Chr Homozygosity Rates)")
    ax[0].set_title("Females (> 0.2 are removed)")
    males.hist(column='F', ax=ax[1])
    ax[1].axvline(0.8, c='r', ls='--')
    ax[1].set_ylabel("Number of Individuals")
    ax[1].set_xlabel("F Value (X Chr Homozygosity Rates)")
    ax[1].set_title("Males (< 0.8 are removed)")
    return fig

def plot_maf_hist(maffile: str="MAF_check.frq"):
    """Plot histograms of minor allele frequency distributions for SNPs.

    Input files should be generated by qc_report.check_sex_report().

    Key arguments:
    --------------
    maffile: str
        prefix for the plink file containing MAF information

    Returns:
    --------
    Figure object
    """
    maf = pd.read_csv(maffile, delimiter=" ", skipinitialspace=True)
    fig = plt.figure(figsize=(8,6))
    plt.hist(maf['MAF'])
    plt.title("MAF Distribution")
    plt.xlabel("MAF")
    plt.ylabel("Number of SNPs")
    return fig

def plot_maf_dropped_hist(maffile: str="MAF_check.frq", threshold: float=0.05):
    """Plot histograms of minor allele frequency (MAF) distributions for SNPs
    for specific frequency thresholds.

    The input file should be generated by the "maf_check_report()" function specified below.

    Key arguments:
    --------------
    maffile: str
        prefix for the plink file containing MAF information
    threshold: float
        MAF threshold

    Returns:
    --------
    Figure object
    """

    maf = pd.read_csv(file, delimiter=" ", skipinitialspace=True)
    rare = maf.loc[maf['MAF'] < threshold]

    fig, ax = plt.subplots(1, 2, figsize=(16,6))
    maf.hist(column='MAF', ax=ax[0])
    ax[0].axvline(threshold, c='r', linestyle='--')
    ax[0].set_title("MAF Distribution of All SNPs\n (< {} are removed)".format(threshold))
    ax[0].set_xlabel("MAF")
    ax[0].set_ylabel("Number of SNPs")
    rare.hist(column='MAF', ax=ax[1])
    ax[1].axvline(threshold, c='r', linestyle='--')
    ax[1].set_title("MAF Distribution of SNPs < {}".format(threshold))
    ax[1].set_xlabel("MAF")
    ax[1].set_ylabel("Number of SNPs")

    # print("dropped SNPs: ", rare.shape[0])
    # print("remaining SNPs: ", maf['SNP'].count() - rare['SNP'].count())
    return fig

def plot_hwe_hist(hwefile: str="plink.hwe", threshold: float=1e-6):
    """Plot histograms of hardy-weinberg equilibrium (HWE) test p-value distributions
    for SNPs.

    The input file should be generated by the qc_report.hardy_weinberg.

    Key arguments:
    --------------
    hwefile: str
        prefix for the plink file containing HWE information
    threshold: float
        p-value threshold

    Returns:
    --------
    Figure object
    """

    hardy = pd.read_csv(file, delimiter=" ", skipinitialspace=True)
    zoomhwe = hardy[hardy['P'] < threshold]

    fig, ax = plt.subplots(1, 2, figsize=(16,6))
    plt.title("HWE Test")
    hardy.hist(column='P', ax=ax[0])
    ax[0].axvline(threshold, c='red', ls='--')
    ax[0].set_xlabel("HWE Exact Test P-Value")
    ax[0].set_ylabel("Number of SNPs")
    ax[0].set_title("HWE P-Value Distribution of All SNPs\n (< {} are removed)".format(threshold))
    zoomhwe.hist(column='P', ax=ax[1])
    ax[1].set_xlabel("HWE Exact Test P-Value")
    ax[1].set_ylabel("Number of SNPs")
    ax[1].set_title("HWE P-Value Distribution of SNPs < {}".format(threshold))
    ax[1].axvline(threshold, c='red', ls='--')
    return fig

def plot_het_hist(het_check_df: pd.DataFrame):
    """Plot histogram of heterozygosity rate distributions for all samples.

    The input file should be generated by the qc_report.heterozygosity_samples().

    Key arguments:
    --------------
    het_check_df: pd.DataFrame
        pandas dataframe containing heterozygosity information

    Returns:
    --------
    Figure object
    """

    fig = plt.figure(figsize=(8,6))
    plt.hist(het_check_df['het_rate'])
    plt.axvline(het_check_df['low_limit'][0], c='red', ls='--')
    plt.axvline(het_check_df['up_limit'][0], c='red', ls='--')
    plt.xlabel("Heterozygosity Rate")
    plt.ylabel("Number of Samples")
    plt.title("Heterozygosity Distribution of All Samples\n (< {:.3f} or > {:.3f} are removed)".format(het_check_df['low_limit'][0], het_check_df['up_limit'][0]))
    return fig

def plot_relatedness_scatter(relatfile: str):
    """Plot histogram of heterozygosity rate distributions for all samples.

    The input file should be generated by the qc_report.relatedness_check().

    Key arguments:
    --------------
    relatfile: str
        file containing relatedness information

    Returns:
    --------
    Figure object
    """
    from matplotlib.ticker import FormatStrFormatter

    relat = pd.read_csv(infile, delimiter=" ", skipinitialspace=True)
    po_z_list = relat.loc[relat['RT'] == "PO"][['Z0', 'Z1']]
    un_z_list = relat.loc[relat['RT'] == "UN"][['Z0', 'Z1']]

    fig = plt.figure(figsize=(8,6))
    if po_z_list.empty and un_z_list.empty:
        return None
    elif po_z_list.empty:
        plt.scatter(x=un_z_list['Z0'], y=un_z_list['Z1'], s=100)
    elif un_z_list.empty:
        plt.scatter(x=po_z_list['Z0'], y=po_z_list['Z1'], s=100)
    else:
        plt.scatter(x=un_z_list['Z0'], y=un_z_list['Z1'], s=100)
        plt.scatter(x=po_z_list['Z0'], y=po_z_list['Z1'], s=100)
    plt.legend(['PO', 'UN'])
    plt.title("Z0 vs Z1 Values for Related (PO) and Unrelated (UN) Individuals")
    return fig

def _make_autopct(values: int):
    """Convert values to percentages.

    Key arguments:
    --------------
    values: int
        values to convert

    Returns:
    --------
    Function object
    """
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct


def write_fail_file(ids_failed: dict, outfile: str="failed_ids"):
    """Write either failed sample IDs or SNPs to two files.

    One is a human-readable csv file, the other is a csv file intended to be
    read into a Pandas DataFrame.

    Key arguments:
    --------------
    ids_failed: dict
        dictionary containing sample and SNPS that failed QC
    outfile: str
        name of file to write failed samples/SNPs IDs to

    Returns:
    --------

    """

    hr = filename + "_hr.csv"
    pr = filename + "_pr.csv"
    with open(hr, "w") as f:
        for k,v in ids_failed.items():
            f.write(str(k) + ": " + str(v) + "\n")
    pd.DataFrame.from_dict(data=ids_failed, orient='index').to_csv(pr, header=True)

def sample_failed_report(imissfile: str="plink.imiss", lmissfile: str="plink.lmiss",
                         sexcheckfile: str="plink.sexcheck",
                         hetfailedfile: str="het_fail_ind.txt",
                         ibdfile: str="pihat_min0.2_in_founders.genome",
                         write: bool=True):
    """Write report of all samples that failed the sample QC checks.

    Key arguments:
    --------------
    imissfile: str
        file containing sample missingness information
    lmissfile: str
        file containing SNP missingness information
    sexcheckfile: str
        file containing sex discrepancy information
    hetfailedfile: str
        file containing heterozygosity information
    ibdfile: str
        file containing cryptic relatedness information (IBD check)
    write: bool
        option to write out ID of samples that failed QC

    Returns:
    --------
    Figure object
    """

    ids = {}
    ids_list = []

    # SNP missingness
    imiss = pd.read_csv(imissfile,  delimiter=" ", skipinitialspace=True)
    lmiss = pd.read_csv(lmissfile,  delimiter=" ", skipinitialspace=True)
    ind_missing_filtered = calculate_missingness(imiss, 'F_MISS', 0.2)
    missing_ids = get_sample_ids(imiss, 'IID', ind_missing_filtered)

    ids['missing'] = missing_ids.tolist()
    ids_list.append(missing_ids.tolist())

    # mismatched sex
    sex = pd.read_csv("plink.sexcheck", delimiter=" ", skipinitialspace=True)
    sex_mismatches = sex.loc[sex['STATUS'] == "PROBLEM"]
    sex_mismatches_counts = sex['STATUS'].value_counts()
    sex_mismatches_ids = sex_mismatches['IID'].tolist()

    ids['sex_mismatches'] = sex_mismatches_ids
    ids_list.append(sex_mismatches_ids)

    # outlying heterozygosity
    het_failed = pd.read_csv("het_fail_ind.txt", delimiter=" ")
    het_failed_ids = het_failed['IID'].tolist()

    ids['het_failed'] = het_failed_ids
    ids_list.append(het_failed_ids)

    # high IBD - pi_hat threshold
    ibd = pd.read_csv("pihat_min0.2_in_founders.genome", delimiter=" ", skipinitialspace=True)
    ibd_ids = ibd['IID1'].tolist()

    ids['relatedness_failed'] = ibd_ids
    ids_list.append(ibd_ids)

    # graph everything
    tests = ['SNP Missingness', 'Sex Mismatches', 'Outlying Heterozygosity', 'Cryptic Relatedness']
    fail_counts = [len(missing_ids), len(sex_mismatches_ids), len(het_failed_ids), len(ibd_ids)]
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

def snps_failed_report(write: bool=False, miss_threshold: float=0.2,
                      maf_threshold: float=0.00001, hwe_threshold: float=1e-6,
                      lmissfile: str="plink.lmiss", maffile: str="MAF_check.frq",
                      hwefile: str="plink.hwe"):
    """Write report of all SNPs that failed QC checks.

    Key arguments:
    --------------
    write: bool
        option to write failed SNPs to file
    miss_threshold: float
        missingness rate threshold for filtering out SNPs
    maf_threshold: float
        MAF threshold for filtering out SNPs
    hwe_threshold: float
        HWE threshold for filtering out SNPs
    lmissfile: str
        file containing missingness information for SNPs
    maffile: str
        file containing MAF information for SNPs
    hwefile: str
        file containing HWE information for SNPs

    Returns:
    --------
    Figure object
    """
    snps = {}
    ids_list = []
    lmiss = pd.read_csv(lmiss_file, delimiter=" ", skipinitialspace=True)

    missing_snps = lmiss.loc[lmiss['F_MISS'] > miss_threshold]
    snps['missing_snps'] = missing_snps['SNP'].tolist()
    ids_list.append(missing_snps['SNP'].tolist())
    # print("total missing snps failed: ", len(missing_snps['SNP'].tolist()))

    # MAF
    maf = pd.read_csv(maf_file, delimiter=" ", skipinitialspace=True)
    rare = maf.loc[maf['MAF'] < maf_threshold]
    snps['maf'] = rare['SNP'].tolist()
    ids_list.append(rare['SNP'].tolist())
    # print("total maf snps failed: ", len(rare['SNP'].tolist()))

    # HWE departures
    hardy = pd.read_csv(hwe_file, delimiter=" ", skipinitialspace=True)
    hwe_failed = hardy.loc[hardy['P'] < hwe_threshold]
    snps['hwe'] = hwe_failed['SNP'].tolist()
    ids_list.append(hwe_failed['SNP'].tolist())
    # print("total hwe snps failed: ", len(hwe_failed['SNP'].tolist()))

    # graph everything
    tests = ['SNP Missingness', 'Minor Allele Frequency', 'Outlying HWE']
    fail_counts = [len(missing_snps['SNP'].tolist()), len(rare['SNP'].tolist()), len(hwe_failed['SNP'].tolist())]
    total_fails = set(x for l in ids_list for x in l)
    # print("total fails: ", len(total_fails))

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

def save_pdf(bfile: str, plots: list):
    """Write a list of matplotlib Figure objects to a pdf file.

    Key arguments:
    --------------
    bfile: str
        prefix for input bfile to perform QC on
    plots: list
        list containing plots generated by the qc_plot module

    Returns:
    --------
    """

    title = infile + "_qc_report.pdf"

    with PdfPages(title) as pdf:
        for plot in plots:
            pdf.savefig(plot)
