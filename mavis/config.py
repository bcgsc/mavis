import argparse
import os
from copy import copy as _copy
from typing import Dict, Optional

import snakemake
import tab
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate as snakemake_validate

from .annotate.file_io import ReferenceFile
from .bam import stats
from .bam.cache import BamCache
from .constants import PROTOCOL, SUBCOMMAND, float_fraction
from .util import WeakMavisNamespace, bash_expands, filepath

CONVERT_OPTIONS = WeakMavisNamespace()
CONVERT_OPTIONS.add(
    'assume_no_untemplated',
    True,
    defn='assume that if not given there is no untemplated sequence between the breakpoints',
)


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


def validate_config(config: Dict, bam_stats: Optional[bool] = False, stage: str = '') -> None:
    """
    Check that the input JSON config conforms to the expected schema as well
    as the other relevant checks such as file exsts
    """
    schema = 'config' if stage != SUBCOMMAND.OVERLAY else 'overlay'

    try:
        snakemake_validate(
            config,
            os.path.join(os.path.dirname(__file__), f'schemas/{schema}.json'),
            set_default=True,
        )
    except Exception as err:
        short_msg = '. '.join(
            [line for line in str(err).split('\n') if line.strip()][:3]
        )  # these can get super long
        raise WorkflowError(short_msg)

    required = []
    if (
        stage not in {SUBCOMMAND.CONVERT}
        or stage == SUBCOMMAND.CLUSTER
        and not config['cluster.uninformative_filter']
    ):
        required.append('reference.annotations')

    if stage == SUBCOMMAND.VALIDATE:
        required.extend(['reference.aligner_reference', 'reference.reference_genome'])

    for req in required:
        if req not in config:
            raise WorkflowError(f'missing required property: {req}')

    if schema == 'config':
        conversion_dir = os.path.join(config['output_dir'], 'converted_outputs')
        # check all assignments are conversions aliases or existing files
        for libname, library in config['libraries'].items():
            assignments = []
            for i, assignment in enumerate(library['assign']):
                if assignment in config.get('convert', {}):
                    # replace the alias with the expected output path
                    converted_output = os.path.join(conversion_dir, f'{assignment}.tab')
                    assignments.append(converted_output)
                elif (
                    not os.path.exists(assignment) and os.path.dirname(assignment) != conversion_dir
                ):
                    raise FileNotFoundError(f'cannot find the expected input file {assignment}')
                else:
                    assignments.append(assignment)
            library['assign'] = assignments

            if not config['skip_stage.validate'] and stage in {
                SUBCOMMAND.VALIDATE,
                SUBCOMMAND.SETUP,
            }:
                if not library.get('bam_file', None) or not os.path.exists(library['bam_file']):
                    raise FileNotFoundError(
                        f'missing bam file for library ({libname}), it is a required input when the validate stage is not skipped'
                    )
                # calculate the bam_stats if the have not been given
                missing_stats = any(
                    [
                        col not in library
                        for col in ['median_fragment_size', 'read_length', 'stdev_fragment_size']
                    ]
                )
                if missing_stats and bam_stats:
                    library.update(calculate_bam_stats(config, libname))

        # expand and check the input files exist for any conversions
        for conversion in config.get('convert', {}).values():
            expanded = []
            for input_file in conversion['inputs']:
                expanded.extend(bash_expands(input_file))
            conversion['inputs'] = expanded

    # make sure all the reference files specified exist and overload with environment variables where applicable
    for ref_type in list(config.keys()):
        if not ref_type.startswith('reference.'):
            continue
        expanded = []
        for input_file in config[ref_type]:
            expanded.extend(bash_expands(input_file))
        config[ref_type] = expanded


def get_metavar(arg_type):
    """
    For a given argument type, returns the string to be used for the metavar argument in add_argument

    Example:
        >>> get_metavar(bool)
        '{True,False}'
    """
    if arg_type in [bool, tab.cast_boolean]:
        return '{True,False}'
    elif arg_type in [float_fraction, float]:
        return 'FLOAT'
    elif arg_type == int:
        return 'INT'
    elif arg_type == filepath:
        return 'FILEPATH'
    return None
