import subprocess
import sys

plink_file="./plink"

def snp_genotype_filter(bfile: str, threshold: float, outfile: str):
    command = "./plink --bfile {} --geno {} --silent --make-bed --out {}".format(bfile, threshold, outfile)
    os.system(command)

def run_plink(bfile, *flags):
    command = f'{plink_file} --bfile {bfile} --silent --make_bed'
    print(flags)
    for flag in flags:
        # command += " {}".format(flag)
        command += f' {flag}'
    print(command)

#%%

run_plink("test", "--geno 0.2", "--make_bed")
