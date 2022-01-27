import argparse
from copy import copy as _copy
from typing import Dict

from .annotate.file_io import ReferenceFile
from .bam import stats
from .bam.cache import BamCache
from .constants import PROTOCOL, float_fraction
from .util import cast_boolean, filepath


def calculate_bam_stats(config: Dict, library_name: str) -> Dict:
    """
    Calculate the read stats for a library from a given bam file
    """
    library = config['libraries'][library_name]
    annotations = ReferenceFile('annotations', *config['reference.annotations'])

    if library['protocol'] == PROTOCOL.TRANS:
        if annotations is None or annotations.is_empty():
            raise AttributeError(
                'missing required attribute: annotations. Annotations must be given for transcriptomes'
            )
        annotations.load()
    bam = BamCache(library['bam_file'], stranded=library['strand_specific'])
    if library['protocol'] == PROTOCOL.TRANS:
        bam_stats = stats.compute_transcriptome_bam_stats(
            bam,
            annotations=annotations.content,
            sample_size=config['bam_stats.sample_size'],
            sample_cap=config['bam_stats.sample_cap'],
            distribution_fraction=config['bam_stats.distribution_fraction'],
        )
        return {
            'median_fragment_size': int(bam_stats.median_fragment_size),
            'read_length': int(bam_stats.read_length),
            'stdev_fragment_size': int(bam_stats.stdev_fragment_size),
            'strand_specific': bam_stats.stranded,
            'strand_determining_read': bam_stats.strand_determining_read,
        }
    bam_stats = stats.compute_genome_bam_stats(
        bam,
        sample_size=config['bam_stats.sample_size'],
        sample_bin_size=config['bam_stats.sample_bin_size'],
        sample_cap=config['bam_stats.sample_cap'],
        distribution_fraction=config['bam_stats.distribution_fraction'],
    )
    return {
        'median_fragment_size': int(bam_stats.median_fragment_size),
        'read_length': int(bam_stats.read_length),
        'stdev_fragment_size': int(bam_stats.stdev_fragment_size),
    }


class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """
    subclass the default help formatter to stop default printing for required arguments
    """

    def _format_args(self, action, default_metavar):
        if action.metavar is None:
            action.metavar = get_metavar(action.type)
        if isinstance(action, RangeAppendAction):
            return '%s' % self._metavar_formatter(action, default_metavar)(1)
        return super(CustomHelpFormatter, self)._format_args(action, default_metavar)

    def _get_help_string(self, action):
        if action.required:
            return action.help
        return super(CustomHelpFormatter, self)._get_help_string(action)

    def add_arguments(self, actions):
        # sort the arguments alphanumerically so they print in the help that way
        actions = sorted(actions, key=lambda x: getattr(x, 'option_strings'))
        super(CustomHelpFormatter, self).add_arguments(actions)


class RangeAppendAction(argparse.Action):
    """
    allows an argument to accept a range of arguments
    """

    def __init__(self, nmin=1, nmax=None, **kwargs):
        kwargs.setdefault('nargs', '+')
        kwargs.setdefault('default', [])
        argparse.Action.__init__(self, **kwargs)
        self.nmin = nmin
        self.nmax = nmax
        assert nmin is not None

    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, [])
        items = _copy(getattr(namespace, self.dest))
        items.append(values)
        if self.nmax is None:
            if len(values) < self.nmin:
                raise argparse.ArgumentError(
                    self, 'must have at least {} arguments. Given: {}'.format(self.nmin, values)
                )
        elif not self.nmin <= len(values) <= self.nmax:
            raise argparse.ArgumentError(
                self, 'requires {}-{} arguments. Given: {}'.format(self.nmin, self.nmax, values)
            )
        setattr(namespace, self.dest, items)


def add_bamstats_to_config(config: Dict):
    """
    Check that the input JSON config conforms to the expected schema as well
    as the other relevant checks such as file exsts
    """
    # check all assignments are conversions aliases or existing files
    for libname, library in config['libraries'].items():
        # calculate the bam_stats if the have not been given
        if any(
            [
                col not in library
                for col in ['median_fragment_size', 'read_length', 'stdev_fragment_size']
            ]
        ):
            library.update(calculate_bam_stats(config, libname))


def get_metavar(arg_type):
    """
    For a given argument type, returns the string to be used for the metavar argument in add_argument

    Example:
        >>> get_metavar(bool)
        '{True,False}'
    """
    if arg_type in [bool, cast_boolean]:
        return '{True,False}'
    elif arg_type in [float_fraction, float]:
        return 'FLOAT'
    elif arg_type == int:
        return 'INT'
    elif arg_type == filepath:
        return 'FILEPATH'
    return None
