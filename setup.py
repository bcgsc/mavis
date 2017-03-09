from setuptools import setup
from mavis import __version__

setup(
    name='MAVIS',
    version=__version__,
    url='https://svn.bcgsc.ca/svn/SVIA/mavis',
    packages=['mavis'],
    scripts=['bin/run_mavis.py', 'bin/mavis_overlay.py'],
    install_requires=[
        'colour',
        'networkx',
        'numpy',
        'biopython',
        'svgwrite',
        'Sphinx',  # for building the documentation only
        'sphinx-rtd-theme',  # for building the documentation only
        'pysam'
    ],
    author_email='creisle@bcgsc.ca',
    dependency_links=[
        'svn+https://svn.bcgsc.ca/svn/SVIA/TSV/tags/v3.1.0#egg=TSV',
        'svn+https://svn.bcgsc.ca/svn/SVIA/vocab/tags/v1.0.0#egg=vocab'
    ],
    test_suite='nose.collector',
    tests_require=['nose']
)
