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
    try:
        import vocab
        if vocab.__version__ != vocab_version:
            raise ImportError('wrong version', vocab.__version__)
        print('Using', vocab.__path__[0])
        print(vocab.__package__ + '==' + vocab.__version__)
        print()
    except ImportError:
        print(['install', '-e', vocab_link, '--trusted-host', '*.bcgsc.ca', '--exists-action', 's'])
        pip.main(['install', '-e', vocab_link, '--trusted-host', '*.bcgsc.ca', '--exists-action', 's'])
    try:
        import TSV
        if TSV.__version__ != tsv_version:
            raise ImportError('wrong version', TSV.__version__)
        print('Using', TSV.__path__[0])
        print(TSV.__package__ + '==' + TSV.__version__)
        print()
    except ImportError:
        print(['install', '-e', tsv_link, '--trusted-host', '*.bcgsc.ca', '--exists-action', 's'])
        pip.main(['install', '-e', tsv_link, '--trusted-host', '*.bcgsc.ca', '--exists-action', 's'])

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
        'pysam>=0.9',
        'TSV=={}'.format(tsv_version),
        'vocab=={}'.format(vocab_version),
        'numpy>=1.11.2',
        'pyvcf==0.6.8',
        'braceexpand==0.1.2'
    ],
    author_email='creisle@bcgsc.ca',
    dependency_links=[
        vocab_link,
        tsv_version
    ],
    setup_requires=['nose>=1.0'],
    test_suite='nose.collector',
    tests_require=['nose', 'timeout-decorator==0.3.3', 'coverage==4.2']
)
