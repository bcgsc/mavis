"""
holds submodules related to structural variants
"""
import subprocess
import os


def get_version():
    v = subprocess.check_output('cd {}; git describe'.format(os.path.dirname(__file__)), shell=True)
    v = v.decode('UTF8')
    v = v.strip()
    return v

__version__ = get_version()
