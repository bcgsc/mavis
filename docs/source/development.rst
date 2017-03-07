Development
===================================

Guidelines for Contributors
-------------------------------

- In general, follow `pep8 <https://www.python.org/dev/peps/pep-0008/>`_ style guides using a maximum line width of 120 characters
- docstrings should follow `sphinx google code style <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_
- any column name which may appear in any of the intermediate or final output files must be defined in :class:`~mavis.constants.COLUMNS`

Formatting Types in docstrings
.................................

if you want to be more explicit with nested types, the following conventions are used throughout the code

- dictionary: ``d = {<key>: <value>}`` becomes ``dict of <value> by <key>``
- list: ``l = [1, 2, 3]`` becomes ``list of int``
- mixed: ``d = {'a': [1, 2, 3], 'b': [4, 5, 6]}`` becomes ``dict of list of int by str``
- tuples: ``('a', 1)`` becomes ``tuple of str and int``


Unit Tests
.................................

- all new code must have unit tests in the tests/ subdirectory
- in general for assertEqual statements, the expected value is given first

Unit tests for this package are run with nose and should be run from the top level directory. You will
need to have installed the python dependencies before running tests. The following command will run
the unit tests and compute the test coverage for the main package.

.. code-block:: bash

    nosetests --with-coverage --cover-html --cover-html-dir=coverage --cover-package=mavis --cover-erase

Virtual Environment
........................

If you are developing in python setting up with a virtual environment can be incredibly helpful. This can be
used to generate the requirements.txt file that pip uses for install. Instructions for setting up the environment
are below

.. code-block:: bash

    pip install virtualenv
    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt

then if new modules/dependencies are required, they can be installed to the virtual python. Before committing
regenerate the requirements.txt file using pip

.. code-block:: bash

    pip freeze > requirements.txt

Major Assumptions
...................

Some assumptions have been made when developing this project. The major ones have been listed here to
facilitate debugging/development if any of these are violated in the future.

- The input bam reads have stored the sequence wrt to the positive/forward strand and have not stored the reverse
  complement.
- The distribution of the fragment sizes in the bam file approximately follows a normal distribution.


Current Limitations
.....................

- Assembling contigs will always fail for repeat sequences as we do not resolve this. Unlike traditional assemblies
  we cannot assume even input coverage as we are taking a select portion of the reads to assemble.
- Currently no attempt is made to group/pair single events into complex events.
- Transcriptome validation uses a collapsed model of all overlapping transcripts and is not isoform specific. Allowing
  for isoform specific validation would be computationally expensive but may be considered as an optional setting for
  future releases.


MAVIS Package Documentation
----------------------------

.. automodule:: mavis
    :special-members: __and__, __or__, __xor__, __len__, __sub__, __add__
    :members:
    :undoc-members:
    :show-inheritance:

.. toctree::
    :maxdepth: -1

    mavis.annotate
    auto/mavis.assemble
    mavis.bam
    auto/mavis.blat
    auto/mavis.breakpoint
    auto/mavis.cluster
    auto/mavis.constants
    auto/mavis.error
    mavis.illustrate
    auto/mavis.interval
    auto/mavis.pairing
    mavis.pipeline
    mavis.validate


Development Goals
-------------------------------

Features to be implemented

.. todolist::

.. |TOOLNAME| replace:: **MAVIS**
