import pandas as pd


def write_txt_file(output_file: str, eids: list):
    with open(output_file, "w") as f:
        for id in eids:
            f.write(str(id) + ",")


def read_txt_file(file: str):
    lines = []
    with open(file, "r") as f:
        for line in f.read().split(','):
            if line != "":
                lines.append(line)
    return lines


def generate_phenofile(ids_file: str, pheno_outfile: str="phenotypes.txt"):
    """Generates phenotype file from text file containing IDs

    Key arguments:
    --------------
    ids_file: str
        file containing list of IDs
    Returns:
    --------
    phenos: list
        list of IDs used to create phenotype file

    """
    ids = read_txt_file(ids_file)
    eids = [id.strip() for id in ids]
    with open(pheno_outfile, "w") as f:
        for id in eids:
            if id:
                line = f"{id} {id} \n"
                f.writelines(line)
    return eids
