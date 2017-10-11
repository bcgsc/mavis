import argparse
from configparser import ConfigParser, ExtendedInterpolation
import os
import re
import warnings

import tab

from . import __version__
from .align import SUPPORTED_ALIGNER
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .bam.cache import BamCache
from .bam.stats import compute_genome_bam_stats, compute_transcriptome_bam_stats
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .constants import DISEASE_STATUS, PROTOCOL
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .tools import SUPPORTED_TOOL
from .util import bash_expands, cast, devnull, ENV_VAR_PREFIX, get_env_variable, MavisNamespace, WeakMavisNamespace
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS


SCHEDULE_DEFAULTS = WeakMavisNamespace(
    queue='transabyss.q',
    validation_memory=16,
    trans_validation_memory=18,
    annotation_memory=12,
    memory=12)


REFERENCE_DEFAULTS = WeakMavisNamespace(
    annotations='',
    reference_genome='',
    template_metadata='',
    masking='',
    aligner_reference='',
    dgv_annotation=''
)


class LibraryConfig(MavisNamespace):

    def __init__(
        self, library, protocol, disease_status, bam_file, inputs, read_length, median_fragment_size,
        stdev_fragment_size, strand_specific, strand_determining_read=2,
        **kwargs
    ):
        self.library = library
        self.protocol = PROTOCOL.enforce(protocol)
        self.bam_file = bam_file
        self.read_length = int(read_length)
        self.median_fragment_size = int(median_fragment_size)
        self.stdev_fragment_size = int(stdev_fragment_size)
        self.strand_specific = cast(strand_specific, bool)
        self.strand_determining_read = int(strand_determining_read)
        self.disease_status = DISEASE_STATUS.enforce(disease_status)
        try:
            self.inputs = [f for f in re.split('[;\s]+', inputs) if f]
        except TypeError:
            self.inputs = inputs

        acceptable = {}
        acceptable.update(CLUSTER_DEFAULTS.__dict__)
        acceptable.update(VALIDATION_DEFAULTS.__dict__)
        acceptable.update(ANNOTATION_DEFAULTS.__dict__)

        for attr, value in kwargs.items():
            setattr(self, attr, cast(value, type(acceptable[attr])))

        if 'MAVIS_FETCH_METHOD_INDIVIDUAL' not in os.environ and 'fetch_method_individual' not in kwargs:
            if self.protocol == PROTOCOL.TRANS:
                self.fetch_method_individual = False
            else:
                self.fetch_method_individual = True

    def flatten(self):
        result = {}
        result.update(self.__dict__)
        result['inputs'] = '\n'.join(result['inputs'])
        return result

    def is_trans(self):
        return True if self.protocol == PROTOCOL.TRANS else False

    @staticmethod
    def build(
        library, protocol, bam_file, inputs,
        annotations=None,
        log=devnull,
        distribution_fraction=0.98,
        sample_cap=3000,
        sample_bin_size=1000,
        sample_size=500,
        **kwargs
    ):
        PROTOCOL.enforce(protocol)

        if protocol == PROTOCOL.TRANS and annotations is None:
            raise AttributeError(
                'missing required attribute: annotations. Annotations must be given for transcriptomes')
        bam = BamCache(bam_file)
        if protocol == PROTOCOL.TRANS:
            bamstats = compute_transcriptome_bam_stats(
                bam,
                annotations=annotations,
                sample_size=sample_size,
                sample_cap=sample_cap,
                distribution_fraction=distribution_fraction
            )
        elif protocol == PROTOCOL.GENOME:
            bamstats = compute_genome_bam_stats(
                bam,
                sample_size=sample_size,
                sample_bin_size=sample_bin_size,
                sample_cap=sample_cap,
                distribution_fraction=distribution_fraction
            )
        else:
            raise ValueError('unrecognized value for protocol', protocol)
        log(bamstats)

        return LibraryConfig(
            library=library, protocol=protocol, bam_file=bam_file, inputs=inputs,
            median_fragment_size=bamstats.median_fragment_size,
            stdev_fragment_size=bamstats.stdev_fragment_size,
            read_length=bamstats.read_length,
            strand_determining_read=bamstats.strand_determining_read,
            **kwargs
        )


def write_config(filename, include_defaults=False, libraries=[], conversions={}, log=devnull):
    """
    Args:
        filename (str): path to the output file
        include_defaults (bool): True if default parameters should be written to the config, False otherwise
        libraries (list of LibraryConfig): library configuration sections
        conversions (dict of list by str): conversion commands by alias name
        log (function): function to pass output logging to
    """
    config = {}

    config['reference'] = REFERENCE_DEFAULTS.flatten()
    for filetype, fname in REFERENCE_DEFAULTS.items():
        if fname is None:
            warnings.warn('filetype {} has not been set. This must be done manually before the configuration file is used'.format(filetype))

    if libraries:
        for lib in libraries:
            config[lib.library] = lib.flatten()

    if include_defaults:
        config['schedule'] = SCHEDULE_DEFAULTS.flatten()
        config['validation'] = VALIDATION_DEFAULTS.flatten()
        config['cluster'] = CLUSTER_DEFAULTS.flatten()
        config['annotation'] = ANNOTATION_DEFAULTS.flatten()
        config['illustrate'] = ILLUSTRATION_DEFAULTS.flatten()
        config['summary'] = SUMMARY_DEFAULTS.flatten()

    config['convert'] = {}
    for alias, command in conversions.items():
        config['convert'][alias] = '\n'.join(command)

    for sec in config:
        for tag, value in config[sec].items():
            if '_regex_' in tag:
                config[sec][tag] = re.sub('\$', '$$', config[sec][tag])
            elif isinstance(value, list):
                config[sec][tag] = '\n'.join([str(v) for v in value])
            else:
                config[sec][tag] = str(value)

    conf = ConfigParser()
    for sec in config:
        conf[sec] = {}
        for tag, val in config[sec].items():
            conf[sec][tag] = val
    log('writing:', filename)
    with open(filename, 'w') as configfile:
        conf.write(configfile)


def validate_and_cast_section(section, defaults, use_defaults=False):
    """
    given a dictionary of values, returns a new dict with the values casted to their appropriate type or set
    to a default if the value was not given
    """
    d = {}
    if use_defaults:
        d.update(defaults.items())
    for attr, value in section.items():
        if attr not in defaults:
            raise KeyError('tag not recognized', attr)
        else:
            d[attr] = cast(value, type(defaults[attr]))
    return d


class MavisConfig:

    def __init__(self, **kwargs):

        # section can be named schedule or qsub to support older versions
        self.schedule = MavisNamespace(**validate_and_cast_section(
            kwargs.pop('schedule', kwargs.pop('qsub', {})), SCHEDULE_DEFAULTS, True
        ))

        # set the global defaults
        for sec, defaults in [
            ('pairing', PAIRING_DEFAULTS),
            ('summary', SUMMARY_DEFAULTS),
            ('validation', VALIDATION_DEFAULTS),
            ('annotation', ANNOTATION_DEFAULTS),
            ('illustrate', ILLUSTRATION_DEFAULTS),
            ('cluster', CLUSTER_DEFAULTS),
            ('reference', REFERENCE_DEFAULTS)
        ]:
            v = MavisNamespace(**validate_and_cast_section(kwargs.pop(sec, {}), defaults, True))
            setattr(self, sec, v)

        SUPPORTED_ALIGNER.enforce(self.validation.aligner)

        for attr, fname in self.reference.items():
            if not os.path.exists(fname):
                raise KeyError(attr, 'file at', fname, 'does not exist')

        # set the conversion section
        self.convert = kwargs.pop('convert', {})
        for attr, val in self.convert.items():
            val = [v for v in re.split('[;\s]+', val) if v]
            if val[0] == 'convert_tool_output':
                if len(val) < 3 or val[2] not in SUPPORTED_TOOL:
                    raise UserWarning(
                        'conversion using the built-in convert_tool_output requires specifying the input file and '
                        'tool name currently supported tools include:', SUPPORTED_TOOL.values())
                inputs = bash_expands(val[1])
                if len(inputs) < 1:
                    raise OSError('input file(s) do not exist', val[1])
                if len(val) == 4:
                    val[3] = tab.cast_boolean(val[3])
                elif len(val) > 4:
                    raise UserWarning(
                        'conversion using the built-in convert_tool_output takes at most 3 arguments')
            self.convert[attr] = val
        self.convert = MavisNamespace(**self.convert)

        # now add the library specific sections
        self.libraries = {}

        for libname, val in kwargs.items():  # all other sections already popped
            d = {}
            d.update(self.cluster.items())
            d.update(self.validation.items())
            d.update(self.annotation.items())
            d.update(val)
            d['library'] = libname
            val['library'] = libname
            self.libraries[libname] = LibraryConfig(**val)
            # now try building the LibraryConfig object
            try:
                lc = LibraryConfig(**d)
                self.libraries[libname] = lc
            except TypeError as terr:  # missing required argument
                try:
                    lc = LibraryConfig.build(**d)
                    self.libraries[libname] = lc
                except Exception as err:
                    raise UserWarning('could not build configuration section for library', libname, err, terr)

    def has_transcriptome(self):
        return any([l.is_trans() for l in self.libraries.values()])

    @staticmethod
    def read(filepath):
        """
        reads the configuration settings from the configuration file

        Args:
            filepath (str): path to the input configuration file

        Returns:
            class:`list` of :class:`Namespace`: namespace arguments for each library
        """
        parser = ConfigParser(interpolation=ExtendedInterpolation())
        parser.read(filepath)
        config_dict = {}

        # get the library sections and add the default settings
        for sec in parser.sections():
            config_dict.setdefault(sec, {}).update(parser[sec].items())
        return MavisConfig(**config_dict)


def add_semi_optional_argument(argname, success_parser, failure_parser, help_msg='', metavar=None):
    """
    for an argument tries to get the argument default from the environment variable
    """
    env_name = ENV_VAR_PREFIX + argname.upper()
    help_msg += ' The default for this argument is configured by setting the environment variable {}'.format(env_name)
    if os.environ.get(env_name, None):
        success_parser.add_argument('--{}'.format(argname), required=False, default=os.environ[env_name], help=help_msg, metavar=metavar)
    else:
        failure_parser.add_argument('--{}'.format(argname), required=True, help=help_msg, metavar=metavar)


def float_fraction(f):
    try:
        f = float(f)
    except ValueError:
        raise argparse.ArgumentTypeError('Argument must be a value between 0 and 1')
    if f < 0 or f > 1:
        raise argparse.ArgumentTypeError('Argument must be a value between 0 and 1')
    return f


def get_metavar(arg_type):
    if arg_type in [bool, tab.cast_boolean]:
        return '{True,False}'
    elif arg_type in [float_fraction, float]:
        return 'FLOAT'
    elif arg_type == int:
        return 'INT'
    return None


def augment_parser(arguments, parser, semi_opt_parser=None, required=None):
    """
    Adds options to the argument parser. Separate function to facilitate the pipeline steps
    all having a similar look/feel
    """
    if required is None:
        required = bool(parser.title.startswith('required'))

    for arg in arguments:

        if arg == 'help':
            parser.add_argument('-h', '--help', action='help', help='show this help message and exit')
        elif arg == 'version':
            parser.add_argument(
                '-v', '--version', action='version', version='%(prog)s version ' + __version__,
                help='Outputs the version number')
        elif arg == 'annotations':
            add_semi_optional_argument(
                arg, semi_opt_parser, parser, 'Path to the reference annotations of genes, transcript, exons, domains, etc.', 'FILEPATH')
        elif arg == 'reference_genome':
            add_semi_optional_argument(arg, semi_opt_parser, parser, 'Path to the human reference genome fasta file.', 'FILEPATH')
        elif arg == 'template_metadata':
            add_semi_optional_argument(arg, semi_opt_parser, parser, 'File containing the cytoband template information.', 'FILEPATH')
        elif arg == 'masking':
            add_semi_optional_argument(arg, semi_opt_parser, parser, metavar='FILEPATH')
        elif arg == 'aligner_reference':
            add_semi_optional_argument(
                arg, semi_opt_parser, parser, 'path to the aligner reference file used for aligning the contig sequences.', 'FILEPATH')
        elif arg == 'dgv_annotation':
            add_semi_optional_argument(
                arg, semi_opt_parser, parser, 'Path to the dgv reference processed to look like the cytoband file.', 'FILEPATH')
        elif arg == 'config':
            parser.add_argument('config', 'path to the config file', metavar='FILEPATH')
        elif arg == 'bam_file':
            parser.add_argument('--bam_file', help='path to the input bam file', required=required, metavar='FILEPATH')
        elif arg == 'read_length':
            parser.add_argument(
                '--read_length', type=int, help='the length of the reads in the bam file',
                required=required, metavar=get_metavar(int))
        elif arg == 'stdev_fragment_size':
            parser.add_argument(
                '--stdev_fragment_size', type=int, help='expected standard deviation in insert sizes',
                required=required, metavar=get_metavar(int))
        elif arg == 'median_fragment_size':
            parser.add_argument(
                '--median_fragment_size', type=int, help='median inset size for pairs in the bam file', required=required,
                metavar=get_metavar(int))
        elif arg == 'library':
            parser.add_argument('--library', help='library name', required=required)
        elif arg == 'protocol':
            parser.add_argument('--protocol', choices=PROTOCOL.values(), help='library protocol', required=required)
        elif arg == 'disease_status':
            parser.add_argument(
                '--disease_status', choices=DISEASE_STATUS.values(), help='library disease status', required=required)
        else:
            value_type = None
            help_msg = None
            default_value = None
            choices = None
            if arg == 'aligner':
                choices = SUPPORTED_ALIGNER.values()
                help_msg = 'aligner to use for aligning contigs'
            if arg == 'strand_specific':
                value_type = tab.cast_boolean
                help_msg = 'indicates that the input is strand specific'
            if arg == 'uninformative_filter':
                help_msg = 'If flag is False then the clusters will not be filtered based on lack of annotation'

            # get default values
            for nspace in [
                    CLUSTER_DEFAULTS,
                    VALIDATION_DEFAULTS,
                    ANNOTATION_DEFAULTS,
                    ILLUSTRATION_DEFAULTS,
                    PAIRING_DEFAULTS,
                    SUMMARY_DEFAULTS]:
                if arg in nspace:
                    default_value = nspace[arg]
                    value_type = type(default_value) if not isinstance(default_value, bool) else tab.cast_boolean
                    if not help_msg:
                        help_msg = 'see user manual for desc'
                    break

            if arg in [
                'assembly_min_remap_coverage',
                'assembly_min_remap_coverage',
                'assembly_strand_concordance',
                'blat_min_identity',
                'contig_aln_min_query_consumption',
                'min_anchor_match',
                'assembly_min_uniq',
                'mask_opacity',
                'min_domain_mapping_match'
            ]:
                value_type = float_fraction

            if help_msg is None:
                raise KeyError('invalid argument', arg)

            parser.add_argument(
                '--{}'.format(arg), choices=choices, metavar=get_metavar(value_type),
                help=help_msg, required=required, default=default_value, type=value_type
            )
