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
    """Generates phenotype file from text file containing IDs that are cases.

    Key arguments:
    --------------
    ids_file: str
        file containing list of IDs that are cases
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


def generate_phenofile_fromfam(ids_file: str, fam_file: str, pheno_outfile: str="phenotypes.txt"):
    """Generates phenotype file from .fam file and file containing a list of
    cases.

    Produces a phenotype file where the first column is FID, second column is IID,
    and third column is 0 if control and 1 if case.

    Key arguments:
    --------------
    ids_file: str
        file containing list of IDs that are cases
    fam_file: str
        .fam file
    Returns:
    --------
    phenos: list
        list of IDs used to create phenotype file

    """
    ids = read_txt_file(ids_file)
    eids = [id.strip() for id in ids]
    fam = pd.read_csv(fam_file, delimiter = " ", usecols = [0, 1], names = ['fid', 'iid'])
    # fam['pheno'] = fam['iid'].apply(lambda x: '1' if x in eids else '0')
    famcopy = fam.copy()
    famcopy['pheno'] = np.where((famcopy['iid'].isin(eids)), 1, 0)
    famcopy.to_csv(pheno_outfile, sep=" ", index=False, header=False)
    return eids
