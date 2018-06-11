import os

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def get_data(*paths):
    return os.path.join(DATA_DIR, *paths)
