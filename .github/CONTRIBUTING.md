
## Getting Started

If you are new to the project a good way to get started is by adding to the documentation, or adding unit tests where
there is a lack of code coverage.

## Install (for Development)

Clone the repository and switch to the development branch

```bash
git clone https://github.com/bcgsc/mavis.git
cd mavis
git checkout develop
```

Set up a python virtual environment. If you are developing in python setting up with a virtual environment can be
incredibly helpful as it allows for a clean install to test. Instructions for setting up the environment
are below

```bash
python3 -m venv venv
source venv/bin/activate
```

Install the MAVIS python package. Running the setup in develop mode will ensure that your code changes are run when you
run MAVIS from within that virtual environment

```bash
pip install -e .[dev]
```

Run the tests and compute code coverage

```bash
pytest tests
```

## Build the Sphinx Documentation

```bash
pip install .[docs]
sphinx-build docs/source/ html
```

The contents of the user manual can then be viewed by opening the build/html/index.html in any available
web browser (i.e. google-chrome, firefox, etc.)


## Deploy to PyPi

Install deployment dependencies

```bash
pip install .[deploy]
```

Build the distribution files

```bash
python setup.py install sdist bdist_wheel
```

Use twine to upload

```bash
twine upload -r pypi dist/*
```


### Reporting a Bug

Please make sure to search through the issues before reporting a bug to ensure there isn't already an open issue.


### Coding Conventions

#### Formatting/Style

- In general, follow [pep8](https://www.python.org/dev/peps/pep-0008/) style guides (except maximum line width)
- docstrings should follow [sphinx google code style](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)
- any column name which may appear in any of the intermediate or final output files must be defined in ``mavis.constants.COLUMNS``


##### Types in docstrings

if you want to be more explicit with nested types, the following conventions are used throughout the code

- dictionary: ``d = {<key>: <value>}`` becomes ``dict of <value> by <key>``
- list: ``l = [1, 2, 3]`` becomes ``list of int``
- mixed: ``d = {'a': [1, 2, 3], 'b': [4, 5, 6]}`` becomes ``dict of list of int by str``
- tuples: ``('a', 1)`` becomes ``tuple of str and int``


#### Tests

- all new code must have unit tests in the tests subdirectory
- in general for assertEqual statements, the expected value is given first

Tests can be run as follows

```bash
pytest tests
```

To run the tests with tox (multiple python installs tested). Note that you will need to have multiple python installs on your path

```bash
tox
```
