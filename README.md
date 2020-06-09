# pygen
python toolkit for running quality control processing and association analysis on genotype call data.

__install using pip for python 3:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pygen.git`

__install from a custom branch:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pygen.git@<branch>`

$ git clone https://github.ibm.com/aur-genomics/pygen.git
$ python setup.py install

### Usage

This module has two submodules:
1. quality_control
2. association

#### Quality control

This submodule is split into two classes:
1. QcSnps
2. QcSamples

Each class expects PLINK binary files as input (bfile - see below). The user should specify the path and name of binary file prefix (e.g HapMap_3_r3_1) as input to the functions, as well as some thresholds where appropriate. Most functions output a figure to analyse the results of the QC step.

In order to import these classes, run the following commands:

$ from pygen import quality_control

Instantiate the classes to run the functions within each class on PLINK data:

$ qc_snps = quality_control.QcSnps()
$ qc_samples = quality_control.QcSamples()

Run QC steps using the available functions:

$ snps_missing_fig = qc_samples.check_snp_missingness(bfile=binary_file_prefix, snp_missingness_cutoff=0.01)

An example QC pipeline script is provided: qc_ukbb.py.

### Testing
An example for a testable function can be found in `exampleTest.py`. This function enables type checking as well as docstring testing.

#### Testing for correctness using the docstring
Each function should include a docstring that contains information about all in- and output variables as well as a test that should pass.

From the command line, navigate to the root directory (`pygen`) and run:
```
pytest --doctest-modules
```

#### Testing typing
From the command line, navigate to the root directory (`pygen`) and run:
```
mypy pygen
```


### Contributing
* Clone the repository
* Checkout a new branch: `git checkout -b my-new-function` (where `my-new-function` is whatever meaningful branch name you come up with)
* Start working on that branch and test, and add, and commit until the new function is ready for review
* If you're creating a new file, don't forget to include it in `pygen/__init__.py`. And if your function has a dependency, include it in `setup.py`.
* Push the branch to the repo: `git push --set-upstream origin my-new-function` (where `my-new-function` is whatever meaningful branch name you came up with)
* From the master branch on github, click on `create Pull Request` and nominate a reviewer
* As a reviewer... review, maybe test. Then comment on the pull request or accept it.
