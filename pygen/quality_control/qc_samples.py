import pandas as pd
import matplotlib.pyplot as plt
from . import plot
from . import report
from . import filter

class QcSamples:
    def __init__(self):
        self.figures = []
        # default names for files
        # miss_out = "plink"
        # sample_out = "sample_missingness_filtered"
        # sexcheck_out = "plink.sexcheck"
        # sex_discrepancy_filtered = "sex_discrepancy_filtered"
        # ld_out = "ld_check"
        # snp_file = "independent_snps"
        # het_filtered = "heterozygosity_filtered"
        # relatedness_filtered = "relatedness_filtered"

    def check_snp_missingness(self, bfile: str, miss_out: str="plink", bfile_out: str="sample_missingness_filtered", snp_missingness_cutoff: float=0.2):
        print("checking snp missingness")
        report.missingness_report(bfile=bfile, outfile=miss_out)
        missing_figs = plot.plot_missingness_hist(miss_file=miss_out)
        self.figures.append(missing_figs)
        print("filter individuals with missingness over threshold {}".format(snp_missingness_cutoff))
        filter.ind_genotype_filter(bfile=bfile, cutoff=snp_missingness_cutoff, outfile=bfile_out)
        return missing_figs

    def check_sex_discrepancy(self, bfile: str="sample_missingness_filtered", sexcheck_out: str="plink.sexcheck", bfile_out: str="sex_discrepancy_filtered"):
        print("check sex discrepancy")
        report.check_sex(bfile=bfile)
        problems_df = report.check_sex_report(file=sexcheck_out)
        check_sex_figs = plot.plot_sexcheck_hist(file=sexcheck_out)
        self.figures.append(check_sex_figs)
        sex_discrepancy = "sex_discrepancy.txt"
        print("filter individuals with sex discrepancies")
        filter.remove_sex(bfile=bfile, remove_file=sex_discrepancy, out=bfile_out)
        return check_sex_figs

    def check_heterozygosity_rate(self, bfile: str="sex_discrepancy_filtered", snp_file: str="independent_snps", ld_out: str="ld_check", window: int=50, shift: int=5, correlation_threshold: float=0.2, correlation_method: str="pairwise", bfile_out: str="heterozygosity_filtered"):
        print("check heterozygosity rates")
        print("get list of non-(highly) correlated snps")
        filter.ld_pruning_filter(bfile=bfile, snp_file=snp_file, outfile=ld_out, window=window, shift=shift, correlation_threshold=correlation_threshold, correlation_method=correlation_method)

        ld_in = ld_out + ".het"
        het_failed = "heterozygosity_failed.txt"
        het_check_df = report.heterozygosity_ind_report(infile=ld_in, outfile=het_failed)
        het_check_fig = plot.plot_het_hist(het_check_df=het_check_df)
        self.figures.append(het_check_fig)
        print("filter individuals with strongly deviating heterozygosity rates")
        hetero_filtered = filter.heterozygosity_filter(bfile=bfile, failed_file=het_failed, outfile=bfile_out)
        return het_check_fig

    def check_cryptic_relatedness(self, bfile: str="heterozygosity_filtered", snp_file: str="independent_snps", threshold: float=0.2, bfile_out: str="relatedness_filtered"):
        print("check cryptic relatedness")
        snp_file_in = snp_file + ".prune.in"
        relatedness_out = "pihat_min{}".format(threshold)
        report.relatedness_check_report(bfile=bfile, indep_snp_file=snp_file_in, outfile=relatedness_out, threshold=threshold)
        relatedness_out_name = relatedness_out + ".genome"
        relat_figs = plot.plot_relatedness_scatter(infile=relatedness_out_name)
        if relat_figs:
            self.figures.append(relat_figs)
            missingness_out = "related_missingness"
            low_call_out = "related_low_call_rate.txt"
            report.missingness_report(bfile=bfile, outfile=missingness_out)
            selected = report.relatives_low_call_rate_report(missingness_out, relatedness_out_name, low_call_out)
            print("filtering related individuals with lowest call rates")
            filter.relatedness_filter(bfile=bfile, remove_file=low_call_out, outfile=bfile_out)
        else:
            filter.rename_filter(bfile=bfile, outfile=bfile_out)
        return relat_figs

    def samples_failed_gen_report(self, bfile: str, write: bool=True, snp_missingness_cutoff: float=0.2, imiss_file: str="plink.imiss", lmiss_file: str="plink.lmiss", sexcheck_file: str="plink.sexcheck", ibd_threshold: float=0.2, ibd_file: str="pihat_min0.2.genome"):
        print("genating sample qc report")
        ibd_file = "pihat_min{}.genome".format(ibd_threshold)
        het_failed_file = "heterozygosity_failed.txt"
        sample_failed_fig = report.sample_failed_report(write=write, miss_threshold=snp_missingness_cutoff, imiss_file=imiss_file, lmiss_file=lmiss_file, sexcheck_file=sexcheck_file, het_failed_file=het_failed_file, ibd_file=ibd_file)
        self.figures.append(sample_failed_fig)
        report_file = bfile + "_samples_qc"
        report.save_pdf(report_file, self.figures)
