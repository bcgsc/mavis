import glob
import os
import shutil

import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


long_running_test = pytest.mark.skipif(
    os.environ.get('RUN_FULL') != '1',
    reason='Only running FAST tests subset',
)

bwa_only = pytest.mark.skipif(not shutil.which('bwa'), reason='missing the command')
blat_only = pytest.mark.skipif(not shutil.which('blat'), reason='missing the command')
todo = pytest.mark.skip(reason='TODO')


def package_relative_file(*paths):
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', *paths))


def get_data(*paths):
    return os.path.join(DATA_DIR, *paths)


def glob_exists(*pos, strict=False, n=1):
    globexpr = os.path.join(*pos)
    file_list = glob.glob(globexpr)
    if strict and len(file_list) == n:
        return file_list[0] if len(file_list) == 1 else file_list
    elif not strict and len(file_list) > 0:
        return file_list
    else:
        print(globexpr)
        print(file_list)
        return False


def glob_not_exists(*pos):
    globexpr = os.path.join(*pos)
    file_list = glob.glob(globexpr)
    return not file_list
