import os
import re

from setuptools import setup


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
