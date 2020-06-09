import pandas as pd
import ukbb_cohort as bb
import os


def filtermainqc(pathmain: str, datafields: list) -> list:
    # Genotype-related QC fields
    # dataFields should include ['22006', '22027', '22021', '22019', '22001', '31']
    main = bb.utils.getColumns(pathToMain=pathmain, dataFields=datafields)
    # only white british ancestry
    causian_main = main[main['22006-0.0'] == '1']
    # exclude heterozygosity and missingness outliers
    het_inliers = causian_main[causian_main['22027-0.0'] != '1']
    # exclude kinship outliers
    kinship_inliers = het_inliers[het_inliers['22021-0.0'] != '-1']
    # exclude relatives
    kinship_filtered = kinship_inliers[kinship_inliers['22021-0.0'] == '0']
    # exclude sex aneuploidies
    sex_filtered = kinship_filtered[kinship_filtered['22019-0.0'] != '1']
    # exclude inviduals with discrepancies between genetic sex and self-reported sex
    sex_dis_filtered = sex_filtered[sex_filtered['31-0.0'] == sex_filtered['22001-0.0']]

    plink_file = sex_dis_filtered[['eid', 'eid']]
    plink_file.to_csv("includeEids.txt", sep=" ", index=False, header=False)
    return sex_dis_filtered['eid']
