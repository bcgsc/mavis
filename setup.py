import os
import re

from setuptools import find_packages, setup

VERSION = '2.2.11'


def parse_md_readme():
    """
    pypi won't render markdown. After conversion to rst it will still not render unless raw directives are removed
    """
    try:
        from m2r import parse_from_file

        rst_lines = parse_from_file('README.md').split('\n')
        long_description = [
            '.. image:: http://mavis.bcgsc.ca/images/acronym.svg\n\n|\n'
        ]  # backup since pip can't handle raw directives
        i = 0
        while i < len(rst_lines):
            if re.match(r'^..\s+raw::.*', rst_lines[i]):
                i += 1
                while re.match(r'^(\s\s+|\t|$).*', rst_lines[i]):
                    i += 1
            else:
                long_description.append(re.sub('>`_ ', '>`__ ', rst_lines[i]))  # anonymous links
                i += 1
        long_description = '\n'.join(long_description)
    except (ImportError, OSError):
        long_description = ''
    return long_description


def check_nonpython_dependencies():
    """
    check that the non-python dependencies have been installed.

    Raises:
        OSError: A dependency is not installed
    """
    import shutil

    aligner = (
        os.environ['MAVIS_ALIGNER']
        if 'MAVIS_ALIGNER' in os.environ and os.environ['MAVIS_ALIGNER']
        else 'blat'
    )
    aligner = re.split(r'\s+', aligner)[0]
    pth = shutil.which(aligner)
    if not pth:
        print('WARNING: Aligner is required. Missing executable: {}'.format(aligner))
    else:
        print('Found: aligner at', pth)


# HSTLIB is a dependency for pysam.
# The cram file libraries fail for some OS versions and mavis does not use cram files so we disable these options
os.environ['HTSLIB_CONFIGURE_OPTIONS'] = '--disable-lzma --disable-bz2 --disable-libcurl'


TEST_REQS = [
    'timeout-decorator>=0.3.3',
    'coverage>=4.2',
    'pycodestyle>=2.3.1',
    'pytest',
    'pytest-cov',
]


DOC_REQS = [
    'mkdocs==1.1.2',
    'markdown_refdocs',
    'mkdocs-material==5.4.0',
    'markdown-include',
    'mkdocs-simple-hooks==0.1.2',
]


INSTALL_REQS = [
    'Distance>=0.1.3',
    'Shapely>=1.6.4.post1',
    'biopython>=1.70, <1.78',
    'braceexpand==0.1.2',
    'colour',
    'networkx==1.11.0',
    'numpy>=1.13.1',
    'pysam>=0.9, <=0.15.2',
    'shortuuid>=0.5.0',
    'svgwrite',
]

DEPLOY_REQS = ['twine', 'm2r', 'wheel']


setup(
    name='mavis',
    version='{}'.format(VERSION),
    url='https://github.com/bcgsc/mavis.git',
    download_url='https://github.com/bcgsc/mavis/archive/v{}.tar.gz'.format(VERSION),
    packages=find_packages(exclude=['tests']),
    description='A Structural Variant Post-Processing Package',
    long_description=parse_md_readme(),
    install_requires=INSTALL_REQS,
    extras_require={
        'docs': DOC_REQS,
        'test': TEST_REQS,
        'dev': ['black', 'flake8'] + DOC_REQS + TEST_REQS + DEPLOY_REQS,
        'deploy': DEPLOY_REQS,
        'tools': ['pyensembl', 'simplejson'],
    },
    tests_require=TEST_REQS,
    setup_requires=['pip>=9.0.0', 'setuptools>=36.0.0,<58'],
    python_requires='>=3.2',
    author='Caralyn Reisle',
    author_email='creisle@bcgsc.ca',
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'mavis = mavis.main:main',
            'calculate_ref_alt_counts = tools.calculate_ref_alt_counts:main',
        ]
    },
    project_urls={'mavis': 'http://mavis.bcgsc.ca'},
)
check_nonpython_dependencies()
