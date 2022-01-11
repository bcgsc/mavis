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

## Build the Documentation

```bash
pip install .[docs]
markdown_refdocs mavis -o docs/package --link
mkdocs build
```

The contents of the user manual can then be viewed by opening the build-docs/index.html in any available web browser 
(i.e. google-chrome, firefox, etc.). Future development to build the Markdown files into HTML and start a development 
server to browse the documentation can be done using: 

```bash
mkdocs serve
```

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

## Reporting a Bug

Please make sure to search through the issues before reporting a bug to ensure there isn't
already an open issue.

## Conventions

### Linting

Use [black](https://github.com/psf/black) with strings off and line length 100

```bash
black src/mavis -S -l 100
```

### Docstrings

docstrings should follow [sphinx google code style](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)

if you want to be more explicit with nested types, please follow the same format
used by [python type annotations](https://docs.python.org/3/library/typing.html)

```text
arg1 (List[str]): a list of strings
```

However using proper type annotations is preferred for new code and then only including the
description of the parameter in the docstring and not its type

```python

def some_function(some_arg: List[str]) -> None:
    """
    Args:
        some_arg: this arg does stuff
    """
```

### Output Columns

any column name which may appear in any of the intermediate or final output files must be defined in `mavis.constants.COLUMNS` as well as added to the [columns glossary](../outputs/columns)

### Tests

- all new code must have unit tests in the tests subdirectory

Tests can be run as follows

```bash
pytest tests
```

### Branching Model

If you are working on a large feature, create a base branch for the feature off develop. Generally
these follow the naming pattern

```bash
git checkout -b integration/issue-<number>-<short-name>
```

If you are working on a smaller feature then simply make a feature branch off develop

```bash
git checkout -b feature/issue-<number>-<short-name>
```

Once ready, a PR should be made to develop and review should be requested from the other developers.

Releases are done by creating a release branch off develop

```bash
git checkout -b release/vX.X.X
```

Updating the version number in setup.py in the release branch, and then making a PR to master.
After the PR has been merged to master a tag/release should be created with the release notes
and a PR to merge master back into develop should be made
