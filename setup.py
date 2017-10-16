from setuptools import setup, find_packages
import pip
import sys
import subprocess
import os
import re

# install local svn dependencies
CWD = os.path.dirname(os.path.abspath(__file__))


def check_nonpython_dependencies():
    import shutil
    pth = shutil.which('samtools')
    if not pth:
        raise OSError('samtools is required')
    else:
        print('FOUND: samtools at', pth)
    pth = shutil.which('blat')
    if not pth:
        print('WARNING: missing blat executable. Blat is the default aligner')
        pth = shutil.which('bwa')
        if pth:
            print('FOUND: bwa at', pth)
            print('please use the MAVIS_ALIGNER environment variable or set the aligner in the config to change the default aligner')
    else:
        print('FOUND: blat at', pth)


def pull_branch_from_git():
    command = 'cd {}; git branch'.format(CWD)
    v = subprocess.check_output(command, shell=True)
    v = v.decode('UTF8')
    for l in v.split('\n'):
        l = l.strip()
        if l.startswith('*'):
            l = re.sub('^\*\s*', '', l)
            return l
    raise OSError('could not parse branch name from git')


def pull_version_from_git():
    command = 'cd {}; git describe --long'.format(CWD)
    v = subprocess.check_output(command, shell=True)
    v = v.decode('UTF8')
    v = v.strip()
    m = re.match('^v?(\d+)\.(\d+)\.(\d+)-(\d+)-g\w+$', v)
    if not m:
        raise OSError('could not parse version number from git', v, '^v?(\d+)\.(\d+)\.(\d+)-\d+-g\d+$')
    commit_number = int(m.group(4))
    branch = pull_branch_from_git()
    flag = 'dev' if branch != 'master' else ''
    if commit_number != 0:
        return '{}.{}.{}.{}{}'.format(m.group(1), m.group(2), m.group(3), flag, commit_number)
    else:
        return '{}.{}.{}'.format(m.group(1), m.group(2), m.group(3))


def write_version_file(version):
    vfile = os.path.join(CWD, 'mavis', 'version.py')
    with open(vfile, 'w') as fh:
        print('writing version to:', vfile, version)
        fh.write('__version__ = \'{}\'\n'.format(version))

try:
    VERSION = pull_version_from_git()
except OSError as err:
    VERSION = None
    while VERSION is None:
        print('failed to auto-detect the version number (requires a git repository)')
        inv = input('please enter the mavis version number that will be used for setup: ')
        if re.match('^\d+\.\d+.\d+(\.(\w\w\w)?\d+)?$', inv):
            VERSION = inv
        else:
            print('error: version is not a valid format. Please follow pep8 versioning i.e. 1.1.1, 1.1.1.dev1, 1.1.1.0, etc.')

print('version:', VERSION)


if any([x in sys.argv for x in ['install', 'develop', 'build']]):
    write_version_file(VERSION)

# HSTLIB is a dependency for pysam
os.environ['HTSLIB_CONFIGURE_OPTIONS'] = '--disable-lzma'  # only required for CRAM files

setup(
    name='mavis',
    version=VERSION,
    url='https://svn.bcgsc.ca/svn/SVIA/mavis',
    packages=find_packages(),
    install_requires=[
        'docutils <0.13.1',
        'colour',
        'networkx==1.11.0',
        'svgwrite',
        'sphinx==1.6.3',  # for building the documentation only
        'sphinx-rtd-theme==0.2.5b1',  # for building the documentation only
        'pysam>=0.9',
        'tab>=0.0.1',
        'numpy>=1.13.1',
        'pyvcf==0.6.8',
        'braceexpand==0.1.2',
        'biopython>=1.70',
        'Distance>=0.1.3'
    ],
    python_requires='>=3',
    author_email='creisle@bcgsc.ca',
    dependency_links=[
        'git+https://svn.bcgsc.ca/bitbucket/scm/prod/tab.git#egg=tab-0.0.1'
    ],
    setup_requires=[
        'nose==1.3.7',
        'numpy>=1.13.1'  # put here b/c biopython doesn't declare this as a setup dependency properly
    ],
    test_suite='nose.collector',
    tests_require=['nose', 'timeout-decorator==0.3.3', 'coverage==4.2', 'mock>=2.0.0'],
    entry_points={'console_scripts': ['mavis = mavis.main:main']}
)
check_nonpython_dependencies()
