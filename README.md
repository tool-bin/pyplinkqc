# pygen
nuts ml/flow functions for genomics

__install using pip for python 3:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pygen.git`

__install from a custom branch:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pygen.git@<branch>`


### Usage

import pygen

... actual functionality coming soon


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
