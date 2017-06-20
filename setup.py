from setuptools import setup
from mavis import __version__
import pip
import sys

tsv_version = '3.1.3'
tsv_link = 'svn+https://svn.bcgsc.ca/svn/SVIA/TSV/tags/v{0}#egg=TSV-{0}'.format(tsv_version)
vocab_version = '1.0.0'
vocab_link = 'svn+https://svn.bcgsc.ca/svn/SVIA/vocab/tags/v{0}#egg=vocab-{0}'.format(vocab_version)

if any([x in sys.argv for x in ['install', 'develop']]):
    # install the svn dependencies. setuptools has tunnel error but pip can do this
    pip.main(['install', '-e', vocab_link])
    pip.main(['install', '-e', tsv_link])

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
        'TSV=={}'.format(tsv_version),
        'vocab=={}'.format(vocab_version),
        'numpy==1.11.2',
        'pyvcf==0.6.8'
    ],
    author_email='creisle@bcgsc.ca',
    dependency_links=[
        vocab_link,
        tsv_version
    ],
    test_suite='nose.collector',
    tests_require=['nose', 'timeout-decorator==0.3.3', 'coverage==4.2']
)
