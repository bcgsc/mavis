from setuptools import setup
from mavis import __version__
import pip
import sys

if any([x in sys.argv for x in ['install', 'develop']]):
    # install the svn dependencies. setuptools has tunnel error but pip can do this
    pip.main(['install', '-e', 'svn+https://svn.bcgsc.ca/svn/SVIA/TSV/tags/v3.1.3#egg=TSV-3.1.3'])
    pip.main(['install', '-e', 'svn+https://svn.bcgsc.ca/svn/SVIA/vocab/tags/v1.0.0#egg=vocab-1.0.0'])

setup(
    name='MAVIS',
    version=__version__,
    url='https://svn.bcgsc.ca/svn/SVIA/mavis',
    packages=['mavis'],
    scripts=['bin/mavis_run.py', 'bin/mavis_overlay.py'],
    install_requires=[
        'docutils <0.13.1',
        'colour',
        'networkx',
        'biopython',
        'svgwrite',
        'Sphinx',  # for building the documentation only
        'sphinx-rtd-theme',  # for building the documentation only
        'pysam==0.9.1.4',
        'TSV==3.1.3',
        'vocab==1.0.0',
        'numpy==1.11.2',
        'pyvcf==0.6.8'
    ],
    author_email='creisle@bcgsc.ca',
    dependency_links=[
        'svn+https://svn.bcgsc.ca/svn/SVIA/TSV/tags/v3.1.3#egg=TSV-3.1.3',
        'svn+https://svn.bcgsc.ca/svn/SVIA/vocab/tags/v1.0.0#egg=vocab-1.0.0'
    ],
    test_suite='nose.collector',
    tests_require=['nose', 'timeout-decorator==0.3.3']
)
