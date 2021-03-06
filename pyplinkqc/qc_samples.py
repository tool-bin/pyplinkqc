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
    qc_report.missingness(bfile, miss_out)
    missing_figs = qc_plot.missingness_hist(miss_out)
    qc_filter.samples_genotypes(bfile, snp_missingness_threshold, bfile_out)
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
    qc_report.run_check_sex(bfile)
    problems_df = qc_report.check_sex(sexcheck_out)
    check_sex_figs = qc_plot.check_sex_hist(sexcheck_out)
    sex_discrepancy = "sex_discrepancy.txt"
    qc_filter.remove_sex(bfile, sex_discrepancy, bfile_out)
    return check_sex_figs

def check_heterozygosity_rate(bfile: str="sex_discrepancy_filtered",
                            snpfile: str="independent_snps", ld_out: str="ld_check",
                            het_out: str="het_check",
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
    qc_filter.ld_pruning(bfile, snpfile, ld_out, window, shift, correlation_threshold,
                      correlation_method)
    qc_report.heterozygosity(ld_out, het_out)
    het_failed = "heterozygosity_failed.txt"
    het_out = het_out + ".het"
    het_check_df = qc_report.heterozygosity_samples(het_out, het_failed)
    het_check_fig = qc_plot.het_hist(het_check_df)
    hetero_filtered = qc_filter.heterozygosity_snps(bfile, het_failed, bfile_out)
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
    relatedness_out = f'pihat_min{threshold}'
    qc_report.relatedness_check(bfile, snpfile_in, relatedness_out, threshold)
    relatedness_out_name = relatedness_out + ".genome"
    relat_figs = qc_plot.relatedness_scatter(relatedness_out_name)
    if relat_figs:
        missingness_out = "related_missingness"
        low_call_out = "related_low_call_rate.txt"
        qc_report.missingness(bfile, missingness_out)
        selected = qc_report.relatives_low_call_rate(missingness_out,
                                                        relatedness_out_name,
                                                        low_call_out)
        qc_filter.relatedness_samples(bfile, low_call_out, bfile_out)
    else:
        qc_filter.rename_filter(bfile, bfile_out)
    return relat_figs

def gen_qc_samples_report(bfile: str, figures_list: list, write: bool=True,
                          snp_missingness_threshold: float=0.2,
                          imissfile: str="plink.imiss",
                          lmissfile: str="plink.lmiss",
                          sexcheckfile: str="plink.sexcheck",
                          ibd_threshold: float=0.2,
                          ibdfile: str="pihat_min0.2.genome"):
    """Generated QC report for samples.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    figures_list: list
        list of figures to write to the report
        (generated by the rest of the functions in this module)
    write: bool
        determines whether to write out failed sample ids
    snp_missingness_threshold: float
        threshold for filtering snps with high missingness rates
    imissfile: str
        file containing samples with high missingness rates
        (generated by check_snps_missingness function)
    lmissfile: str
        file containing snps with high missingness rates
        (generated by check_snp_missingness function)
    sexcheckfile: str
        file containing samples that have sex discrepancies
        (generated by check_sex_discrepancy function)
    ibd_threshold: float
        threshold for removing samples with high ibd
    ibdfile: str
        file containing ibc coefficients for samples
        (generated by check_cryptic_relatedness function)

    Returns:
    --------

    """
    het_failed_file = "heterozygosity_failed.txt"
    sample_failed_fig = qc_report.samples_failed(write, snp_missingness_threshold,
                                                    imissfile, lmissfile, sexcheckfile,
                                                    het_failed_file, ibdfile)
    report_file = bfile + "_samples_qc"
    qc_report.save_pdf(report_file, figures_list)
