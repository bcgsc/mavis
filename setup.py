import os
from setuptools import setup, find_packages
import re


VERSION = '1.6.8'


def parse_md_readme():
    """
    pypi won't render markdown. After conversion to rst it will still not render unless raw directives are removed
    """
    try:
        from m2r import parse_from_file
        rst_lines = parse_from_file('README.md').split('\n')
        long_description = ['.. image:: http://mavis.bcgsc.ca/docs/latest/_static/acronym.svg\n\n']  # backup since pip can't handle raw directives
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
    pth = shutil.which('samtools')
    if not pth:
        raise OSError('Samtools is required. Missing executable: samtools')
    print('Found: samtools at', pth)
    aligner = os.environ['MAVIS_ALIGNER'] if 'MAVIS_ALIGNER' in os.environ and os.environ['MAVIS_ALIGNER'] else 'blat'
    aligner = re.split(r'\s+', aligner)[0]
    pth = shutil.which(aligner)
    if not pth:
        raise OSError('Aligner is required. Missing executable: {}'.format(aligner))
    print('Found: aligner at', pth)

# HSTLIB is a dependency for pysam
os.environ['HTSLIB_CONFIGURE_OPTIONS'] = '--disable-lzma'  # only required for CRAM files


setup(
    name='mavis',
    version='{}'.format(VERSION),
    url='https://github.com/bcgsc/mavis.git',
    download_url='https://github.com/bcgsc/mavis/archive/v{}.tar.gz'.format(VERSION),
    packages=find_packages(exclude=['tests']),
    description='A Structural Variant Post-Processing Package',
    long_description=parse_md_readme(),
    install_requires=[
        'docutils==0.14',
        'colour',
        'networkx==1.11.0',
        'svgwrite',
        'sphinx==1.6.3',  # for building the documentation only
        'sphinx-rtd-theme==0.2.5b1',  # for building the documentation only
        'pysam>=0.9',
        'numpy>=1.13.1',
        'pyvcf==0.6.8',
        'braceexpand==0.1.2',
        'biopython>=1.70',
        'Distance>=0.1.3',
        'setuptools>=36.6.0',
        'shortuuid>=0.5.0',
        'm2r>=0.1.12',
        'Shapely==1.6.4.post1'
    ],
    python_requires='>=3',
    author='Caralyn Reisle',
    author_email='creisle@bcgsc.ca',
    setup_requires=[
        'numpy>=1.13.1',  # put here b/c biopython doesn't declare this as a setup dependency properly
        'setuptools>=36.6.0',
        'nose',
        'timeout-decorator==0.3.3',
        'coverage==4.2',
        'mock>=2.0.0',
        'nose-capturestderr==1.2',
        'nose-exclude>=0.5.0'
    ],
    test_suite='nose.collector',
    entry_points={'console_scripts': ['mavis = mavis.main:main']},
    project_urls={'mavis': 'http://mavis.bcgsc.ca'}
)
check_nonpython_dependencies()
