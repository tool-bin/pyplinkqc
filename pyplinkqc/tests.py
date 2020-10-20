# import qc_samples
from pyplinkqc import qc_samples, qc_snps
import configparser
import subprocess
# import qc_snps
# import association

# This script implements common quality control (QC) on both SNP and sample levels
# and genome-wide association analysis on the QC SNPs and samples.
# The QC results are reported in PDF files found within the "examples" directory

# Specify PLINK binary files prefix
bfile="../examples/HapMap_3_r3_1"
config_file="../plink_config.conf"


# Specify thresholds
snp_missingness_cutoff=0.2
get_autosomal=False
maf_threshold=0.01
hwe_threshold=1e-10

# Perform sample QC procedures
#snps_missing_fig = qc_samples.check_snp_missingness(bfile=bfile, snp_missingness_threshold=snp_missingness_cutoff)

sex_discrepancy_fig = qc_samples.check_sex_discrepancy(bfile="sample_missingness_filtered")

het_rate_fig = qc_samples.check_heterozygosity_rate(bfile="sex_discrepancy_filtered")

related_fig = qc_samples.check_cryptic_relatedness(bfile="heterozygosity_filtered")
figs = [sex_discrepancy_fig, het_rate_fig, related_fig]
# Generate sample QC report
qc_samples.gen_qc_samples_report(bfile=bfile, figures_list=figs, snp_missingness_threshold=snp_missingness_cutoff)


#%%
# Perform SNP QC procedures
from pyplinkqc import qc_snps
# Specify PLINK binary files prefix
bfile="../examples/HapMap_3_r3_1"
config_file="../plink_config.conf"


# Specify thresholds
snp_missingness_cutoff=0.2
get_autosomal=False
maf_threshold=0.01
hwe_threshold=1e-10
snps_missing_fig = qc_snps.check_snp_missingness(bfile=bfile, snp_missingness_threshold=snp_missingness_cutoff)

maf_check, maf_drop = qc_snps.check_maf(bfile="snp_missingness_filtered", get_autosomal=get_autosomal, maf_threshold=maf_threshold)

hwe_check = qc_snps.check_hwe(bfile="maf_filtered", hwe_threshold=hwe_threshold)

figures_list = [snps_missing_fig, maf_check, maf_drop, hwe_check]
# Generate SNP QC report
qc_snps.gen_qc_snps_report(bfile=bfile, figures_list=figures_list, snp_missingness_threshold=snp_missingness_cutoff, maf_threshold=maf_threshold, hwe_threshold=hwe_threshold)

# Perform genome-wide association analysis using logistic regression
#association.perform_cov_assoc(bfile="plink.hwe", outfile="log_association.assoc.log")
