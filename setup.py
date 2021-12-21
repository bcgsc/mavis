import os
import re

from setuptools import setup

# HSTLIB is a dependency for pysam.
# The cram file libraries fail for some OS versions and mavis does not use cram files so we disable these options
os.environ['HTSLIB_CONFIGURE_OPTIONS'] = '--disable-lzma --disable-bz2 --disable-libcurl'

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

setup()
check_nonpython_dependencies()