import numpy as np
import pandas as pd
import argparse
import qc_samples as qcsamples
import qc_snps as qcsnps
import filter
import report
import plot

## This script implements quality control procedures on UKBB genotype data
## It is specific for UKBB genotype data, and follows on from the QC steps performed by the UKBB team
## First we filter the .fam files for individuals that passed the QC procedures perfored by the UKBB team
## With these updated plink files, we perform the standard QC steps as recommended by Ben Neale's group (http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas)

# Guidelines:
#			i. Filter individuals based on UKBB QC procedures
# 				1) Population substructure i.e only white british genetic ancestry
# 				2) Remove closely related individuals (or at lest one of a pair of related individauls)
# 				3) Individuals with sex chromosome aneuploidies
# 				4) Individuals who withdrawn consent
# 			ii. SNP QC:
# 				1) MAF > 0.1%
# 				2) HWE p-value > 1e-10
#               3) INFO score > 0.8 (directly from UKBB)

def qc_all(bfile):
    ## Set parameters to those used in Ben Neale's documentation
    snp_missingness_cutoff=0.2
    get_autosomal=False
    maf_threshold=0.01
    hwe_threshold=1e-10

    qc_samples = qcsamples.QcSamples()
    qc_snps = qcsnps.QcSnps()
    # Only check missingness of SNPs as samples with outliers for missingness have already been filtered. Also the data is stored per chromosome so the proportion of missing SNPs per individuals is not possible to calculate without merging these data. Merging leads to memory issues.
    snps_missing_fig = qc_snps.check_snp_missingness(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

    sex_discrepancy_fig = qc_samples.check_sex_discrepancy(bfile="snp_missingness_filtered")

    maf_check, maf_drop = qc_snps.check_maf(bfile="sex_discrepancy_filtered", get_autosomal=get_autosomal, maf_threshold=maf_threshold)

    hwe_check = qc_snps.check_hwe(bfile="maf_filtered", hwe_threshold=hwe_threshold)

    het_rate_fig = qc_samples.check_heterozygosity_rate(bfile="hwe_filtered")

    related_fig = qc_samples.check_cryptic_relatedness(bfile="heterozygosity_filtered")

    qc_samples.samples_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

    qc_snps.snps_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff, maf_threshold=maf_threshold, hwe_threshold=hwe_threshold)

def qc_samples(bfile):
    ## Set parameters to those used in Ben Neale's documentation
    snp_missingness_cutoff=0.2

    qc_samples = qcsamples.QcSamples()

    snps_missing_fig = qc_samples.check_snp_missingness(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

    sex_discrepancy_fig = qc_samples.check_sex_discrepancy(bfile="snp_missingness_filtered")

    het_rate_fig = qc_samples.check_heterozygosity_rate(bfile="sex_discrepancy_filtered")

    related_fig = qc_samples.check_cryptic_relatedness(bfile="heterozygosity_filtered")

    qc_samples.samples_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

def qc_snps(bfile):
    ## Set parameters to those used in Ben Neale's documentation
    snp_missingness_cutoff=0.2
    get_autosomal=False
    maf_threshold=0.01
    hwe_threshold=1e-10

    qc_snps = qcsnps.QcSnps()

    snps_missing_fig = qc_snps.check_snp_missingness(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

    maf_check, maf_drop = qc_snps.check_maf(bfile="snp_missingness_filtered", get_autosomal=get_autosomal, maf_threshold=maf_threshold)

    hwe_check = qc_snps.check_hwe(bfile="maf_filtered", hwe_threshold=hwe_threshold)

    qc_snps.snps_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff, maf_threshold=maf_threshold, hwe_threshold=hwe_threshold)


def main(bfile: str, filter_file: str="", full: bool=True, sample: bool=False, snp: bool=False):
    if filter_file != "":
        print("filtering plinks files using {}".format(filter_file))
        bfile_out = bfile + "_eid_filtered"
        filter.individual_filter(bfile=bfile, keepfile=filter_file, outfile=bfile_out)
        bfile = bfile_out
    if sample:
        qc_samples(bfile)
    elif snp:
        qc_snps(bfile)
    else:
        qc_all(bfile)


if __name__=="__main__":
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("bfile", help="provide the prefix name for the .bed, .fam files (e.g HapMap_3_r3_1(.bed, .bim, .fam))")
    parser.add_argument("-filter_file", help="option to filter plink files for a specific set of individuals. provide the include file with EIDs as input to this argument", default="")
    #parser.add_argument("-full", help="run full QC analysis on samples and SNPs", action="store_true", default=False)
    parser.add_argument("-only_samples", help="option to only run QC on the samples (individuals)", action="store_true", default=False)
    parser.add_argument("-only_snps", help="option to only run QC on the SNPs", action="store_true", default=False)

    args = parser.parse_args()
    bfile = args.bfile
    filter_file = args.filter_file
    #full_var = args.full
    sample_var = args.only_samples
    snp_var = args.only_snps

    main(bfile, filter_file=filter_file, sample=sample_var, snp=snp_var)
