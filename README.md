# pyplinkqc
Python package for performing quality control (QC) steps and association analysis on genotype call data. The package contains modules that function as wrappers around the **PLINK** (v1.9) command-line tool. In order to use this package, you will need to download the PLINK executable. Please visit the following website, and download the correct version for your operating system:

https://www.cog-genomics.org/plink/1.9/
https://zzz.bwh.harvard.edu/plink/download.shtml

### Installation

pyplinkqc can be installed by running the following commands in the terminal:

<!-- __install using pip for python 3:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pyplinkqc.git`

__install from a specific branch:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pyplinkqc.git@<branch>` -->

__install by cloning branch and running setup.py:__

`$ git clone https://github.ibm.com/aur-genomics/pygen.git`

`$ python setup.py install`

__NB__: We highly recommend the use of a virtual environment when installing the pyplinkqc package. Conda is a great place to start: https://docs.conda.io/projects/conda/en/latest/user-guide/index.html.

### Usage

Once pylinkqc and PLINK have been installed, you will need to configure the path to the PLINK executable. This can be done by updating the "plink.conf" file in the parent directory of this package. The path is specified under the "[PATHS]" section of the configuration file:

`[PATHS]
plink_path="plink"
`

You can update the value of "plink_path" to the appropriate path to your PLINK executable. Once this is done, the package will automatically parse the configuration file, and set the correct path to the executable.

You can now import pyplinkqc within existing python scripts and start using it:

```
import pyplinkqc
bfile_path = "examples/HapMap_3_r3_1"
snp_missing_fig = qc_snps.check_snp_missingness(bfile_path)
```

As you can see, the package assumes the use of PLINK binary files to store genetic data. Examples of such files can be found in the "examples" directory of this repository.

###  Modules

This package has several modules:

1. qc_snps.py - functions to perform QC steps at SNP-level
2. qc_samples.py - functions to perform QC steps at individuals/sample-level
3. qc_filter.py - functions to filter genotype call data
4. qc_report.py - functions to report on QC results
5. qc_plot.py - functions to create plots of QC results
6. association.py - function to perform genome-wide association studies

The following sections outline the intended usage of each of the modules, along with examples for how to run the functions using the example dataset provided in the "examples" directory. Feel free to follow along using your favourite IDE.

#### Quality control

Modules 1 to 5 listed above can be used to perform common quality control procedures on genotype call data, as well as report the results. Quality control procedures can be performed on on a SNP-level, an indvidua/sample-level, or both.

Examples of running quality control steps at SNP-level are given below.

```
from pylinkqc import qc_snps

bfile_path = "examples/HapMap_3_r3_1"

snp_missingness_cutoff=0.2
get_autosomal=False
maf_threshold=0.01
hwe_threshold=1e-10

snp_missing_fig = qc_snps.check_snp_missingness(bfile_path)

maf_check, maf_drop = qc_snps.check_maf(bfile="snp_missingness_filtered", get_autosomal=get_autosomal, maf_threshold=maf_threshold)

hwe_check = qc_snps.check_hwe(bfile="maf_filtered", hwe_threshold=hwe_threshold)

qc_snps.snps_failed_gen_report(bfile=bfile, snp_missingness_cutoff=snp_missingness_cutoff, maf_threshold=maf_threshold, hwe_threshold=hwe_threshold)
```

As shown from the code snippet above, each function that performs a QC step expects the path and name of PLINK binary file prefix (e.g HapMap_3_r3_1).

Screenshots of the generated QC report are shown below:

An example script that implements a full QC pipeline at both SNP and sample level is given in the "examples" directory.

#### Assocation

This module implements basic association testing for GWAS studies. The module implements the following association tests:

1. chi-squared allelic association tests (see [here](https://zzz.bwh.harvard.edu/plink/anal.shtml#cc) for further information)
2. linear or logistic regession tests (see [here](https://zzz.bwh.harvard.edu/plink/anal.shtml#glm) for further information)

An example of performing a GWAS using a logistic regression test is shown below
```
from pyplinkqc import association

bfile_path = "examples/HapMap_3_r3_1"

association.perform_cov_assoc(bfile=bfile_path, outfile="log_association", type="log")
```
The above snippet will run the association analysis on the provided PLINK binary files and save the results in the "log_association.assoc.log" file.

An example script that implements a full QC and association pipeline is given the "examples" directory.

### Contributing

As a collaborator, please create a branch and create a pull request when ready. To contribute otherwise, please fork directory and create pull requests. Github issues are also welcome.

### Citation

If you've found this tool useful in your work, please use the following citations:

pyplinkqc citation:

**pyplinkqc: a Python package for genetic studies**
Nathalie Willems, Isabell Kiral, Benjamin Goudey

PLINK citation:

Package:     PLINK (1.9)
Author:      Shaun Purcell
URL:         http://pngu.mgh.harvard.edu/purcell/plink/

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
PLINK: a toolset for whole-genome association and population-based
linkage analysis. American Journal of Human Genetics, 81.
