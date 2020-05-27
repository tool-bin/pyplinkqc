import pandas as pd
import os

def generate_pheno_plink(bfile: str, phenofile: str, outfile: str):
    """Generates plink binary files annotated with case-control phenotypes.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    phenofile: str
        file that contains the IDs of the individuals that have the phenotype - expects FID is the first column and IID is the second column
    outfile: str
        prefix for the output plink binary files
    """
    command = f"./plink --bfile {bfile} --make-pheno {phenofile} '*' --out {outfile} --make-bed"
    os.system(command)


def perform_simple_assoc(bfile: str, adjust: bool=False, outfile: str):
    """Performs 1 df chi-squared allelic assocation test. Assumes the plink binary files contain phenotype annotations (6th column of the .fam/.ped file). PLINK automatically infers qualitative vs quantitative (contains values other than 0, 1, 2 or missing) phenotypes.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    adjust: bool
        optional parameter specifying whether to employ Bonferroni Correction for adjustment of the p-values
    outfile: str
        prefix for the output plink binary files
    """
    if adjust:
        command = "./plink --bfile {} --adjust --assoc --out {}".format(bfile, outfile)
    else:
        command = "./plink --bfile {} --assoc --out {}".format(bfile, outfile)
    os.system(command)


def perform_cov_assoc(bfile: str, type: str="linear", cov: str="", outfile: str):
    """Performs linear/logistic regression association analysis.

    Key arguments:
    --------------
    bfile: str
        prefix for plink binary files (.bed, .bim, .fam)
    type: str
        specifies whether to run a "linear" or "logistic" regression
    cov: str
        optional parameter that specifies a covariate file (e.g txt file listing a covariate)
    outfile: str
        prefix for the output plink binary files
    """
    if type=="linear":
        command = "./plink --bfile {} --linear --out {}".format(bfile, outfile)
    elif type=="log":
        command = "./plink --bfile {} --log --out {}".format(bfile, outfile)
    else:
        raise Exception("{} is not a supported regression model. Please try linear or log".format(type))
    if cov != "":
        command = command + " --cov {}".format(cov)
    os.system(command)
