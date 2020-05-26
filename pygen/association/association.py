import pandas as pd
import os
import pickle

def read_bfile(binaryfile: str):
    with open(binaryfile, 'rb') as file:
        tmplist = pickle.load(file)
    filelist = [str(i) for i in tmplist]
    return filelist


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



def perform_simple_assoc(bfile: str, phenofile: str, outfile: str):
    command = f"./plink --bfile {bfile} --make-pheno {phenofile} '*' --assoc --out {outfile}"
    os.system(command)

def perform_adjust_assoc(bfile: str, phenofile: str, outfile: str):
    command = f"./plink --bfile {bfile} --make-pheno {phenofile} '*' --assoc --adjust --out {outfile}"
    os.system(command)


#plink --bfile mydata --make-pheno p1.list * --assoc

#%%

bfile="HapMap_3_r3_1"
eidfile = "../../../modellingScripts/nat/cohortPipeline/merged_gpclincal_maindataset_eids_borderline glaucoma.txt"
caselist = generate_phenofile(bfile=bfile, eidfile=eidfile, outfile="caselist.csv")

#%%

perform_assoc(bfile=bfile, phenofile="caselist.csv", outfile="test_assoc")


/Users/nathaliewillems/Box/Projects/Genomics/data/UKBB_Data/main_data/ukb41268_head100.csv
