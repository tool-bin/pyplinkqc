import pandas as pd
import matplotlib.pyplot as plt
from . import plot
from . import report
from . import filter

# default names for files:
# miss_out = "plink"
# snp_filtered = "snp_missingness_filtered"
# maf_check = "MAF_check.frq"
# maf_filtered = "maf_filtered"
# hwe_check = "plink.hwe"
# hwe_filtered = "hwe_filtered"

def check_snp_missingness(bfile: str, miss_out: str="plink", snp_missingness_cutoff: float=0.2, bfile_out: str="snp_missingness_filtered"):
    report.missingness_report(bfile=bfile, outfile=miss_out)
    missing_figs = plot.plot_missingness_hist(miss_file=miss_out)
    # self.figures.append(missing_figs)
    filter.snp_genotype_filter(bfile=bfile, threshold=snp_missingness_cutoff, outfile=bfile_out)
    return missing_figs

def check_maf(bfile: str="snp_missingness_filtered", get_autosomal: bool=False, maf_check: str="MAF_check.frq", maf_threshold: float=0.01, bfile_out: str="maf_filtered"):
    if get_autosomal:
        auto_out = "snp_1_22.txt"
        auto_df = filter.autosomal_snps_filter(bfile=bfile, outfile=auto_out)
        bfile_tmp = "maf_auto"
        filter.select_autosomal_snp(bfile=bfile, auto_file=auto_out, outfile=bfile_tmp)
        bfile = bfile_tmp
    report.maf_check_report(bfile=bfile)
    maf_check_figs = plot.plot_maf_hist(file=maf_check)
    # self.figures.append(maf_check_figs)
    maf_filtered = filter.maf_filter(bfile=bfile, threshold=maf_threshold, outfile=bfile_out)
    maf_drop_figs = plot.plot_maf_dropped_hist()
    # self.figures.append(maf_drop_figs)
    return maf_check_figs, maf_drop_figs

def check_hwe(bfile: str="maf_filtered", hwe_check: str="plink.hwe", hwe_threshold: float=1e-6, control: bool=True, bfile_out: str="hwe_filtered"):
    report.hardy_weinberg_report(bfile=bfile)
    hwe_figs = plot.plot_hwe_hist(file=hwe_check, threshold=hwe_threshold)
    # self.figures.append(hwe_figs)
    filter.hardy_weinberg_filter(bfile=bfile, threshold=hwe_threshold, control=control, outfile=bfile_out)
    return hwe_figs

def snps_failed_gen_report(bfile: str, figures_list: list, write: bool=False, snp_missingness_cutoff: float=0.2, maf_threshold: float=0.01, hwe_threshold: float=1e-6, lmiss_file: str="plink.lmiss", maf_file: str="MAF_check.frq", hwe_file: str="plink.hwe"):
    snps_failed_fig = report.snps_failed_report(miss_threshold=snp_missingness_cutoff, maf_threshold=maf_threshold, hwe_threshold=hwe_threshold, lmiss_file=lmiss_file, maf_file=maf_file, hwe_file=hwe_file)
    # self.figures.append(snps_failed_fig)
    report_file = bfile + "_snps_qc"
    report.save_pdf(report_file, figures_list)
