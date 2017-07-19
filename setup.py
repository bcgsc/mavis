from setuptools import setup, find_packages
import pip
import sys
import subprocess
import os
import re


# install local svn dependencies
cwd = os.path.dirname(os.path.abspath(__file__))
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


def check_nonpython_dependencies():
    from mavis.validate.constants import DEFAULTS
    from mavis.config import get_env_variable
    import shutil
    from mavis.blat import get_blat_version
    from mavis.bam.read import get_samtools_version

    aligner = get_env_variable('aligner', DEFAULTS.aligner, str)
    exe = shutil.which(aligner)
    if not exe:
        raise UserWarning('missing dependency', aligner)
    else:
        print('\nUsing', exe)
        if aligner == 'blat':
            version = get_blat_version()
            print('Current version:', aligner, version)

    tool = 'samtools'
    exe = shutil.which(tool)
    if not exe:
        raise UserWarning('missing dependency', aligner)
    else:
        print('\nUsing', exe)
        print('Current version:', tool, 'v{}.{}.{}'.format(*get_samtools_version()))


def pull_version_from_git():
    command = 'cd {}; git describe --long --all'.format(cwd)
    v = subprocess.check_output(command, shell=True)
    v = v.decode('UTF8')
    v = v.strip()
    m = re.match('^v?(\d+)\.(\d+)\.(\d+)-(\d+)-g\w+$', v)
    if not m:
        raise OSError('could not parse version number from git', v, '^v?(\d+)\.(\d+)\.(\d+)-\d+-g\d+$')
    return '{}.{}.{}'.format(m.group(1), m.group(2), int(m.group(3)) + int(m.group(4)))


version = pull_version_from_git()
print('version:', version)
vfile = os.path.join(cwd, 'mavis', 'version.py')

if any([x in sys.argv for x in ['install', 'develop']]):
    with open(vfile, 'w') as fh:
        print('writing version to:', vfile, version)
        fh.write('__version__ = \'{}\'\n'.format(version))

setup(
    name='mavis',
    version=version,
    url='https://svn.bcgsc.ca/svn/SVIA/mavis',
    packages=find_packages(),
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
    tests_require=['nose', 'timeout-decorator==0.3.3', 'coverage==4.2'],
    entry_points={'console_scripts': ['mavis = mavis.main:main']}
)

check_nonpython_dependencies()
