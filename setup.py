from setuptools import setup, find_packages
import pip
import sys
import subprocess
import os
import re

# install local svn dependencies
cwd = os.path.dirname(os.path.abspath(__file__))
vocab_version = '1.0.0'
vocab_link = 'svn+https://svn.bcgsc.ca/svn/SVIA/vocab/tags/v{0}#egg=vocab-{0}'.format(vocab_version)


def install_vocab():
    try:
        import vocab
        if vocab.__version__ != vocab_version:
            raise ImportError('wrong version', vocab.__version__)
        print('Using', vocab.__path__[0])
        print(vocab.__package__ + '==' + vocab.__version__)
        print()
    except ImportError:
        pip.main(['install', '-e', vocab_link, '--trusted-host', '*.bcgsc.ca', '--exists-action', 's'])


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
    command = 'cd {}; git branch'.format(cwd)
    v = subprocess.check_output(command, shell=True)
    v = v.decode('UTF8')
    for l in v.split('\n'):
        l = l.strip()
        if l.startswith('*'):
            l = re.sub('^\*\s*', '', l)
            return l
    raise OSError('could not parse branch name from git')


def pull_version_from_git():
    command = 'cd {}; git describe --long'.format(cwd)
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
    vfile = os.path.join(cwd, 'mavis', 'version.py')
    with open(vfile, 'w') as fh:
        print('writing version to:', vfile, version)
        fh.write('__version__ = \'{}\'\n'.format(version))

try:
    version = pull_version_from_git()
except OSError as err:
    version = None
    while version is None:
        print('failed to auto-detect the version number (requires a git repository)')
        inv = input('please enter the mavis version number that will be used for setup: ')
        if re.match('^\d+\.\d+.\d+(\.(\w\w\w)?\d+)?$', inv):
            version = inv
        else:
            print('error: version is not a valid format. Please follow pep8 versioning i.e. 1.1.1, 1.1.1.dev1, 1.1.1.0, etc.')

print('version:', version)


if any([x in sys.argv for x in ['install', 'develop', 'build']]):
    install_TSV()
    install_vocab()
    write_version_file(version)

# HSTLIB is a dependency for pysam
os.environ['HTSLIB_CONFIGURE_OPTIONS'] = '--disable-lzma'  # only required for CRAM files

setup(
    name='mavis',
    version=version,
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
        'vocab=={}'.format(vocab_version),
        'numpy>=1.13.1',
        'pyvcf==0.6.8',
        'braceexpand==0.1.2',
        'biopython>=1.70',
        'Distance>=0.1.3'
    ],
    python_requires='>=3',
    author_email='creisle@bcgsc.ca',
    dependency_links=[
        vocab_link,
        'git+https://creisle@svn.bcgsc.ca/bitbucket/scm/prod/tab.git#egg=tab-0.0.1'
    ],
    setup_requires=[
        'nose>=1.0',
        'numpy>=1.13.1'  # put here b/c biopython doesn't declare this as a setup dependency properly
    ],
    test_suite='nose.collector',
    tests_require=['nose', 'timeout-decorator==0.3.3', 'coverage==4.2'],
    entry_points={'console_scripts': ['mavis = mavis.main:main']}
)
check_nonpython_dependencies()
