# nutsgen
nuts ml/flow functions for genomics

__install using pip for python 3:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pygen.git`

__install from a custom branch:__

`pip3 install git+ssh://git@github.ibm.com/aur-genomics/pygen.git@<branch>`


### Usage

import pygen

... actual functionality coming soon


### Testing
An example for a testable function can be found in `test.py`. This function enables type checking as well as docstring testing.

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
