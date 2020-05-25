import pandas as pd
import os
import pickle

def read_bfile(binaryfile: str):
    with open(binaryfile, 'rb') as file:
        tmplist = pickle.load(file)
    filelist = [str(i) for i in tmplist]
    return filelist


def generate_phenofile(bfile: str, eidfile: str, outfile: str):
    caselist = read_bfile(eidfile)
    pheno = pd.DataFrame(columns=['FID', 'IID'])
    pheno['FID'] = caselist
    pheno['IID'] = caselist
    pheno.to_csv(outfile, sep=" ", header=False, index=False)
    return caselist


def perform_assoc(bfile: str, phenofile: str, outfile: str):
    command = "plink --bfile {} --make-pheno {} '*' --assoc --out {}".format(bfile, phenofile, outfile)
    os.system(command)


#plink --bfile mydata --make-pheno p1.list * --assoc

#%%

bfile="HapMap_3_r3_1"
eidfile = "../cohortPipeline/merged_gpclincal_maindataset_eids_borderline glaucoma.txt"
caselist = generate_phenofile(bfile=bfile, eidfile=eidfile, outfile="caselist.csv")

#%%

perform_assoc(bfile=bfile, phenofile="caselist.csv", outfile="test_assoc")


/Users/nathaliewillems/Box/Projects/Genomics/data/UKBB_Data/main_data/ukb41268_head100.csv
