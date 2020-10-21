import pandas as pd
import numpy as np
import os
import subprocess
from .run_plink import run_plink


def snp_genotypes(bfile: str, threshold: float, outfile: str):
    """Filters SNPs based on missing genotype rate.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    threshold: float
        threshold to use to filter SNP genotype
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    run_plink(bfile, f'--geno {threshold}', f'--out {outfile}')

def samples_genotypes(bfile: str, threshold: float, outfile: str):
    """Filters samples based on missing genotype rate.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    threshold: float
        threshold to use to filter individuals based on missing SNPs
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # command = "./plink --bfile {} --mind {} --silent --make-bed --out {}".format(bfile, cutoff, outfile)
    # os.system(command)
    run_plink(bfile, f'--mind {threshold}', f'--out {outfile}')

def individuals(bfile: str, keepfile: str, outfile: str):
    """Filters out individuals based on sample ID

    Filters out individuals based on a space/tab-delimited text file with family IDs in first column and within-family IDs in the second column. Individuals who are not in the file are removed.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    keepfile: str
        path to space/tab-delimited text file
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # command = "./plink --bfile {} --keep {} --silent --make-bed --out {}".format(bfile, keepfile, outfile)
    # os.system(command)
    run_plink(bfile, f'-- keep {keepfile}', f'--out {outfile}')


def impute_sex(bfile: str, outfile: str):
    """Change sex assignment based on imputed values.

    Imputed values are calculated from X chromosome inbreeding coefficients.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # command = "./plink --bfile {} --impute-sex --silent --make-bed --out {}".format(bfile, outfile)
    # os.system(command)
    run_plink(bfile, '--impute sex', f'--out {outfile}')

def remove_sex(bfile: str, removefile: str, outfile: str):
    """Remove individuals with sex discrepancies.

    Expects that removefile was generated by "check_sex" function in qc_report module.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    removefile: str
        file generated by check_sex function in qc_report module
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    if os.path.isfile(removefile):
        # command = "./plink --bfile {} --remove {} --make-bed --out {}".format(bfile, remove_file, out)
        # os.system(command)
        run_plink(bfile, f'--remove {removefile}', f'--out {outfile}')
    else:
        print("error in remove_sex function! {} file is not found!").format(removefile)

def _autosomal_snps_file(bfile: str, outfile: str):
    """Generate file of autosomal SNPs.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    bim_file = bfile + ".bim"
    bim = pd.read_csv(bim_file, delimiter="\t", skipinitialspace=True)
    bim.columns = ['chrom', 'snp', 'cm', 'pos', 'a0', 'a1']
    bim['chrom'] = bim['chrom'].astype('int32')
    bim_filtered = bim.loc[(bim['chrom'] >= 1) & (bim['chrom'] <= 22)]
    bim_filtered['snp'].to_csv(outfile, index=False, sep=' ')
    # return bim_filtered

def autosomal_snp(bfile: str, autofile: str, outfile: str):
    """Filter autosomal SNPs.

    SNPs that are not autosomal are filtered out. Expects file generated by _autosomal_snps_file function (autofile arg).

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    autofile: str
        path to file containg autosomal SNPs (generated by _autosomal_snps_file function)
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # # command = "./plink --bfile {} --extract {} --silent --make-bed --out {}".format(bfile, auto_file, outfile)
    # os.system(command)
    run_plink(bfile, f'--extract {autofile}')

def maf(bfile: str, threshold: float, outfile: str):
    """Filter variants with minor allele frequency below threshold.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    threshold: float
        threshold for minor allele frequency
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # command = "./plink --bfile {} --maf {} --make-bed --out {}".format(bfile, threshold, outfile)
    # os.system(command)
    run_plink(bfile, f'--maf {threshold}', f'--out {outfile}')

def hardy_weinberg_test(bfile: str, control: bool, threshold: float, outfile: str):
    """Filter out variants with HWE exact test p-value below threshold.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    control: bool
        determines whether to consider only controls (True) or controls and cases (False)
    threshold: float
        threshold for hardy-weinberg equilibrium exact test (p-value)
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    if control:
        # command = "./plink --bfile {} --hwe {} --silent --make-bed --out {}".format(bfile, threshold, outfile)
        run_plink(bfile, f'--hwe {threshold}', f'--out {outfile}')
    else:
        # command = "./plink --bfile {} --hwe include-nonctrl {} --silent --make-bed --out {}".format(bfile, threshold, outfile)
        run_plink(bfile, f'--hwe include-nonctrl {threshold}', f'--out {outfile}')
    # os.system(command)
    #os.system(command2)

def ld_pruning(bfile: str, snpfile: str, outfile: str, window: int=50, shift: int=5, correlation_threshold: int=0.2, correlation_method: str="pairwise"):
    """Filter out SNPs in high linkage disequilibirium.

    Only keep SNPs that are in approximate linkage equilibrium with each other.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    snpfile: str
        file to write independent SNPs to
    outfile: str
        prefix to write plink binary files to
    window: int
        number of SNPs considered in a window for the calculation
    shift: int
        number of SNPs to shift the window by
    correlation_threshold: float
        threshold for filter SNPs above this correlation
    correlation_method: str
        method to use for calculating the correlation (default: pairwise)

    Returns:
    --------

    """
    correlation_methods = ['multiple', 'pairwise']
    if correlation_method not in correlation_methods:
        raise ValueError(f'{correlation_method} not a valid choice, please choose from {correlation_methods}')
    if correlation_method=="multiple":
        if threshold < 1:
            raise ValueError(f'{threshold} must be above 1')
        # command = "./plink --bfile {} --indep {} {} {} --out {}".format(bfile, window, shift, correlation_threshold, snp_file)
        run_plink(bfile, f'--indep {window} {shift} {correlation_threshold}', f'--out {snpfile}', make_bed=False)
    elif correlation_method=="pairwise":
        # command = "./plink --bfile {} --indep-pairwise {} {} {} --out {}".format(bfile, window, shift, correlation_threshold, snp_file)
        run_plink(bfile, f'--indep-pairwise {window} {shift} {correlation_threshold}', f'--out {snpfile}', make_bed=False)
    snp_in = snpfile + ".prune.in"
    # command2 = "./plink --bfile {} --extract {} --het --out {}".format(bfile, snp_in, outfile)
    run_plink(bfile, f'--extract {snp_in}', f'--out {outfile}')
    # os.system(command)
    # os.system(command2)

def heterozygosity_snps(bfile: str, failedfile: str, outfile: str):
    """Filter out SNPs with high heterozygosity rates.

    SNPs are removed based on SNPs listed in the failedfile.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    failedfile: str
        file contained list of SNPs to filter out
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # command = "./plink --bfile {} --remove {} --make-bed --out {}".format(bfile, failed_file, outfile)
    # os.system(command)
    run_plink(bfile, f'--remove {failedfile}', f'--out {outfile}')

def relatedness_samples(bfile: str, removefile: str, outfile: str):
    """Filter samples that are related.

    Samples are removed based on removefile.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    removefile: str
        file containing list of samples to filter out
    outfile: str
        prefix for the output plink binary files

    Returns:
    --------

    """
    # command = "./plink --bfile {} --remove {} --make-bed --out {}".format(bfile, remove_file, outfile)
    # os.system(command)
    run_plink(bfile, f'--remove {removefile}', f'--out {outfile}')
