from pyplinkqc import qc_samples
from pyplinkqc import qc_snps
from pyplinkqc import association

# This script implements common quality control (QC) on both SNP and sample levels
# and genome-wide association analysis on the QC SNPs and samples.
# The QC results are reported in PDF files found within the "examples" directory

# Specify PLINK binary files prefix
bfile="HapMap_3_r3_1"

# Specify thresholds
snp_missingness_cutoff=0.2
get_autosomal=False
maf_threshold=0.01
hwe_threshold=1e-10

# Perform sample QC procedures
snps_missing_fig = qc_samples.check_snp_missingness(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

sex_discrepancy_fig = qc_samples.check_sex_discrepancy(bfile="snp_missingness_filtered")

het_rate_fig = qc_samples.check_heterozygosity_rate(bfile="sex_discrepancy_filtered")

related_fig = qc_samples.check_cryptic_relatedness(bfile="heterozygosity_filtered")

# Generate sample QC report
qc_samples.samples_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

# Perform SNP QC procedures
snps_missing_fig = qc_snps.check_snp_missingness(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff)

maf_check, maf_drop = qc_snps.check_maf(bfile="snp_missingness_filtered", get_autosomal=get_autosomal, maf_threshold=maf_threshold)

hwe_check = qc_snps.check_hwe(bfile="maf_filtered", hwe_threshold=hwe_threshold)

# Generate SNP QC report
qc_snps.snps_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff, maf_threshold=maf_threshold, hwe_threshold=hwe_threshold)

# Perform genome-wide association analysis using logistic regression
association.perform_cov_assoc(bfile="plink.hwe", outfile="log_association.assoc.log")
