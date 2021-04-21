import collections
import os

from snakemake.utils import validate as snakemake_validate


class ImmutableDict(collections.Mapping):
    def __init__(self, data):
        self._data = data

    def __getitem__(self, key):
        return self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)


def get_by_prefix(config, prefix):
    return {k.replace(prefix, ''): v for k, v in config.items() if k.startswith(prefix)}


DEFAULTS = {}
snakemake_validate(
    DEFAULTS,
    os.path.join(os.path.dirname(__file__), 'config.json'),
    set_default=True,
)
DEFAULTS = ImmutableDict(DEFAULTS)
