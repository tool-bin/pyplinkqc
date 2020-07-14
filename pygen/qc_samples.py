import pandas as pd
import matplotlib.pyplot as plt
from . import qc_plot
from . import qc_report
from . import qc_filter

# default names for files
# miss_out = "plink"
# sample_out = "sample_missingness_filtered"
# sexcheck_out = "plink.sexcheck"
# sex_discrepancy_filtered = "sex_discrepancy_filtered"
# ld_out = "ld_check"
# snp_file = "independent_snps"
# het_filtered = "heterozygosity_filtered"
# relatedness_filtered = "relatedness_filtered"

def check_snp_missingness(bfile: str, miss_out: str="plink",
                          bfile_out: str="sample_missingness_filtered",
                          snp_missingness_threshold: float=0.2):
    """Filters SNPs with high missingness rates.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    miss_out: str
        file to write missingness report to
    bfile_out: str
        prefix for the output plink binary files
    snp_missingness_threshold: float
        threshold to use for SNPs missingness rate

    Returns:
    --------
        Figure object
    """
    report.missingness_report(bfile, miss_out)
    missing_figs = plot.plot_missingness_hist(miss_out)
    filter.samples_genotypes(bfile, snp_missingness_threshold, bfile_out)
    return missing_figs

def check_sex_discrepancy(bfile: str="sample_missingness_filtered",
                         sexcheck_out: str="plink.sexcheck",
                         bfile_out: str="sex_discrepancy_filtered"):
    """Filters out samples with sex discrepancies.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    sexcheck_out: str
        file to write sexcheck report to
    bfile_out: str
        prefix for the output plink binary files

    Returns:
    --------
        Figure object
    """
    report.check_sex(bfile)
    problems_df = report.check_sex_report(sexcheck_out)
    check_sex_figs = plot.plot_sexcheck_hist(sexcheck_out)
    sex_discrepancy = "sex_discrepancy.txt"
    filter.remove_sex(bfile, sex_discrepancy, bfile_out)
    return check_sex_figs

def check_heterozygosity_rate(bfile: str="sex_discrepancy_filtered",
                            snpfile: str="independent_snps", ld_out: str="ld_check",
                            window: int=50, shift: int=5, correlation_threshold: float=0.2,
                            correlation_method: str="pairwise",
                            bfile_out: str="heterozygosity_filtered"):
    """Filters samples with high heterozygosity rates.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    snpfile: str
        file to write independent SNPs to
    ld_out: str
        file to write linkage disequilibirium results to
    window: int
        number of SNPs considered in a window for the calculation
    shift: int
        number of SNPs to shift the window by
    correlation_threshold: float
        threshold for filter SNPs above this correlation
    correlation_method: str
        method to use for calculating the correlation (default: pairwise)
    bfile_out: str
        prefix for the output plink binary files

    Returns:
    --------
        Figure object
    """
    filter.ld_pruning(bfile, snpfile, ld_out, window, shift, correlation_threshold,
                      correlation_method)
    ld_in = ld_out + ".het"
    het_failed = "heterozygosity_failed.txt"
    het_check_df = report.heterozygosity_ind_report(ld_in, het_failed)
    het_check_fig = plot.plot_het_hist(het_check_df)
    hetero_filtered = filter.heterozygosity_snps(bfile, het_failed, bfile_out)
    return het_check_fig

def check_cryptic_relatedness(bfile: str="heterozygosity_filtered",
                              snpfile: str="independent_snps", threshold: float=0.2,
                              bfile_out: str="relatedness_filtered"):
    """Filter samples with cryptic relatedness.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    snpfile: str
        file containing list of independent SNPs
        (SNPs in linkage equilibrium)
    threshold: float
        pi_hat threshold
    bfile_out: str
        prefix for the output plink binary files

    Returns:
    --------
        Figure object
    """
    snpfile_in = snpfile + ".prune.in"
    # relatedness_out = "pihat_min{}".format(threshold)
    relatedness_out = f'pihat_min{threshold}'
    report.relatedness_check_report(bfile, snp_file_in, relatedness_out, threshold)
    relatedness_out_name = relatedness_out + ".genome"
    relat_figs = plot.plot_relatedness_scatter(relatedness_out_name)
    if relat_figs:
        missingness_out = "related_missingness"
        low_call_out = "related_low_call_rate.txt"
        report.missingness_report(bfile, missingness_out)
        selected = report.relatives_low_call_rate_report(missingness_out,
                                                        relatedness_out_name,
                                                        low_call_out)
        filter.relatedness_samples(bfile, low_call_out, bfile_out)
    else:
        filter.rename_filter(bfile, bfile_out)
    return relat_figs

def gen_qc_samples_report(bfile: str, write: bool=True, snp_missingness_threshold: float=0.2,
imissfile: str="plink.imiss", lmissfile: str="plink.lmiss", sexcheckfile: str="plink.sexcheck",
ibd_threshold: float=0.2, ibdfile: str="pihat_min0.2.genome", figures_list: list=figures):
    """Generated QC report for samples.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    write: bool
        determines whether to write out failed sample IDs
    snp_missingness_threshold: float
        threshold for filtering SNPs with high missingness rates
    imissfile: str
        file containing samples with high missingness rates
        (generated by check_snps_missingness function)
    lmissfile: str
        file containing SNPs with high missingness rates
        (generated by check_snp_missingness function)
    sexcheckfile: str
        file containing samples that have sex discrepancies
        (generated by check_sex_discrepancy function)
    ibd_threshold: float
        threshold for removing samples with high IBD
    ibdfile: str
        file containing IBC coefficients for samples
        (generated by check_cryptic_relatedness function)
    figures_list: list
        list of figures to write to the report
        (generated by the rest of the functions in this module)

    Returns:
    --------

    """
    het_failed_file = "heterozygosity_failed.txt"
    sample_failed_fig = report.sample_failed_report(write, snp_missingness_threshold,
                                                    imissfile, lmissfile, sexcheck_file,
                                                    het_failed_file, ibd_file)
    report_file = bfile + "_samples_qc"
    report.save_pdf(report_file, figures)
