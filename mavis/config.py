import argparse
from configparser import ConfigParser, ExtendedInterpolation
from copy import copy as _copy
import os
import re
import warnings

import tab

from . import __version__
from .align import SUPPORTED_ALIGNER
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .annotate.file_io import load_annotations
from .bam.cache import BamCache
from .bam.stats import compute_genome_bam_stats, compute_transcriptome_bam_stats
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .constants import DISEASE_STATUS, SUBCOMMAND, PROTOCOL, float_fraction
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .submit import OPTIONS as SUBMIT_OPTIONS
from .submit import SCHEDULER
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .tools import SUPPORTED_TOOL
from .util import bash_expands, cast, devnull, ENV_VAR_PREFIX, MavisNamespace, WeakMavisNamespace, get_env_variable, log_arguments
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS


REFERENCE_DEFAULTS = WeakMavisNamespace(
    annotations='',
    reference_genome='',
    template_metadata='',
    masking='',
    aligner_reference='',
    dgv_annotation=''
)

CONVERT_OPTIONS = WeakMavisNamespace()
CONVERT_OPTIONS.add('assume_no_untemplated', True, defn='assume that if not given there is no untemplated sequence between the breakpoints')


class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """
    subclass the default help formatter to stop default printing for required arguments
    """
    def _format_args(self, action, default_metavar):
        if isinstance(action, RangeAppendAction):
            return '%s' % self._metavar_formatter(action, default_metavar)(1)
        return super(CustomHelpFormatter, self)._format_args(action, default_metavar)

    def _get_help_string(self, action):
        if action.required:
            return action.help
        return super(CustomHelpFormatter, self)._get_help_string(action)


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
                    self, 'must have at least {} arguments. Given: {}'.format(self.nmin, values))
        elif not self.nmin <= len(values) <= self.nmax:
            raise argparse.ArgumentError(
                self, 'requires {}-{} arguments. Given: {}'.format(self.nmin, self.nmax, values))
        setattr(namespace, self.dest, items)


def cast_if_not_none(value, cast_type):
    """
    cast a value to a given type unless it is None

    Example:
        >>> cast_if_not_none('1', int)
        1
        >>> cast_if_not_none(None, int)
        None
        >>> cast_if_not_none('null', int)
        None
        >>> cast_if_not_none('', int)
        None
        >>> cast_if_not_none('none', int)
        None
    """
    if value is None or str(value).lower() in ['none', 'null', '']:
        return None
    return cast(value, cast_type)


class LibraryConfig(MavisNamespace):
    """
    holds library specific configuration information
    """
    def __init__(
        self, library, protocol, disease_status, bam_file=None, inputs=None, read_length=None, median_fragment_size=None,
        stdev_fragment_size=None, strand_specific=False, strand_determining_read=2,
        **kwargs
    ):
        MavisNamespace.__init__(self)
        self.library = library
        self.protocol = PROTOCOL.enforce(protocol)
        self.bam_file = bam_file
        self.read_length = cast_if_not_none(read_length, int)
        self.median_fragment_size = cast_if_not_none(median_fragment_size, int)
        self.stdev_fragment_size = cast_if_not_none(stdev_fragment_size, int)
        self.strand_specific = cast(strand_specific, bool)
        self.strand_determining_read = int(strand_determining_read)
        self.disease_status = DISEASE_STATUS.enforce(disease_status)
        try:
            self.inputs = [f for f in re.split(r'[;\s]+', inputs) if f]
        except TypeError:
            self.inputs = inputs if inputs is not None else []

        for attr, value in kwargs.items():
            for namespace in [CLUSTER_DEFAULTS, VALIDATION_DEFAULTS, ANNOTATION_DEFAULTS]:
                if attr not in namespace:
                    continue
                setattr(self, attr, namespace.type(attr)(value))
                break

    def flatten(self):
        result = MavisNamespace.flatten(self)
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
        """
        Builds a library config section and gathers the bam stats
        """
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

    @classmethod
    def parse_args(cls, *args):
        # '<name>', '(genome|transcriptome)', '<diseased|normal>', '[strand_specific]', '[/path/to/bam/file]'
        if len(args) < 4:
            return LibraryConfig(args[0], protocol=args[1], disease_status=args[2])
        elif len(args) < 5:
            return LibraryConfig(args[0], protocol=args[1], disease_status=args[2], strand_specific=args[3])
        return LibraryConfig(args[0], protocol=args[1], disease_status=args[2], strand_specific=args[3], bam_file=args[4])


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
        config['schedule'] = SUBMIT_OPTIONS.flatten()
        config['validate'] = VALIDATION_DEFAULTS.flatten()
        config['cluster'] = CLUSTER_DEFAULTS.flatten()
        config['annotate'] = ANNOTATION_DEFAULTS.flatten()
        config['illustrate'] = ILLUSTRATION_DEFAULTS.flatten()
        config['summary'] = SUMMARY_DEFAULTS.flatten()

    config['convert'] = CONVERT_OPTIONS.flatten()
    for alias, command in conversions.items():
        if alias in CONVERT_OPTIONS:
            raise UserWarning('error in writing config. Alias for conversion product cannot be a setting', alias, CONVERT_OPTIONS.keys())
        config['convert'][alias] = '\n'.join(command)

    for sec in config:
        for tag, value in config[sec].items():
            if '_regex_' in tag:
                config[sec][tag] = re.sub(r'\$', '$$', config[sec][tag])
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


def validate_section(section, namespace, use_defaults=False):
    """
    given a dictionary of values, returns a new dict with the values casted to their appropriate type or set
    to a default if the value was not given
    """
    new_namespace = MavisNamespace()
    if use_defaults:
        for attr, value in namespace.items():
            new_namespace.add(attr, value, cast_type=namespace.type(attr))

    for attr, value in section.items():
        if attr not in namespace:
            raise KeyError('tag not recognized', attr)
        else:
            try:
                new_namespace.add(attr, namespace.type(attr)(value), cast_type=namespace.type(attr))
            except ValueError:
                raise ValueError('failed casting {} with value {}'.format(attr, value))
    return new_namespace


class MavisConfig:

    def __init__(self, **kwargs):

        # section can be named schedule or qsub to support older versions
        self.schedule = validate_section(kwargs.pop('schedule', kwargs.pop('qsub', {})), SUBMIT_OPTIONS, True)

        # set the global defaults
        for sec, defaults in [
            ('pairing', PAIRING_DEFAULTS),
            ('summary', SUMMARY_DEFAULTS),
            ('validate', VALIDATION_DEFAULTS),
            ('annotate', ANNOTATION_DEFAULTS),
            ('illustrate', ILLUSTRATION_DEFAULTS),
            ('cluster', CLUSTER_DEFAULTS),
            ('reference', REFERENCE_DEFAULTS)
        ]:
            setattr(self, sec, validate_section(kwargs.pop(sec, {}), defaults, True))

        SUPPORTED_ALIGNER.enforce(self.validate.aligner)

        for attr, fname in self.reference.items():
            if not os.path.exists(fname):
                raise OSError(attr, 'file at', fname, 'does not exist')

        # set the conversion section
        self.convert = kwargs.pop('convert', {})
        for attr, val in self.convert.items():
            if attr in CONVERT_OPTIONS:
                self.convert[attr] = CONVERT_OPTIONS.type(attr)(val)
                continue
            val = [v for v in re.split(r'[;\s]+', val) if v]
            if not val:
                raise UserWarning('conversion tag requires arguments', attr)
            if val[0] == 'convert_tool_output':
                try:
                    val[-1] = tab.cast_boolean(val[-1])
                except TypeError:
                    val.append(False)
                if len(val) < 4 or val[-2] not in SUPPORTED_TOOL.values():
                    raise UserWarning(
                        'conversion using the built-in convert_tool_output requires specifying the input file(s) and '
                        'tool name. currently supported tools include:', SUPPORTED_TOOL.values(), 'given', val)
                expanded_inputs = []
                for file_expr in val[1:-2]:
                    expanded = bash_expands(file_expr)
                    if not expanded:
                        raise OSError('input file(s) do not exist', val[1:-2])
                    expanded_inputs.extend(expanded)
                val = [val[0]] + expanded_inputs + val[-2:]
            self.convert[attr] = val
        self.convert = MavisNamespace(**self.convert)

        # now add the library specific sections
        self.libraries = {}

        for libname, val in kwargs.items():  # all other sections already popped
            d = {}
            d.update(self.cluster.items())
            d.update(self.validate.items())
            d.update(self.annotate.items())
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
        if not os.path.exists(filepath):
            raise OSError('File does not exist: {}'.format(filepath))
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
        required = required = bool(success_parser.title.startswith('required'))
        success_parser.add_argument('--{}'.format(argname), required=required, default=os.environ[env_name], help=help_msg, metavar=metavar)
    else:
        required = required = bool(failure_parser.title.startswith('required'))
        failure_parser.add_argument('--{}'.format(argname), required=required, help=help_msg, metavar=metavar)


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
    return None


def augment_parser(arguments, parser, semi_opt_parser=None, required=None):
    """
    Adds options to the argument parser. Separate function to facilitate the pipeline steps
    all having a similar look/feel
    """
    if required is None:
        try:
            required = bool(parser.title.startswith('required'))
        except AttributeError:
            pass

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
        elif arg == 'skip_stage':
            parser.add_argument(
                '--skip_stage', choices=[SUBCOMMAND.CLUSTER, SUBCOMMAND.VALIDATE], action='append', default=[],
                help='Use flag once per stage to skip. Can skip clustering or validation or both')
        elif arg == 'strand_specific':
            parser.add_argument(
                '--strand_specific', type=tab.cast_boolean, metavar=get_metavar(bool),
                default=False, help='indicates that the input is strand specific')
        else:
            value_type = None
            help_msg = None
            default_value = None
            choices = None
            if arg == 'aligner':
                choices = SUPPORTED_ALIGNER.values()
                help_msg = 'aligner to use for aligning contigs'
            if arg == 'uninformative_filter':
                help_msg = 'If flag is False then the clusters will not be filtered based on lack of annotation'
            if arg == 'scheduler':
                choices = SCHEDULER.keys()

            # get default values
            for nspace in [
                    CLUSTER_DEFAULTS,
                    VALIDATION_DEFAULTS,
                    ANNOTATION_DEFAULTS,
                    ILLUSTRATION_DEFAULTS,
                    PAIRING_DEFAULTS,
                    SUMMARY_DEFAULTS,
                    SUBMIT_OPTIONS,
                    CONVERT_OPTIONS]:
                if arg in nspace:
                    default_value = nspace[arg]
                    value_type = type(default_value) if not isinstance(default_value, bool) else tab.cast_boolean
                    if not help_msg:
                        help_msg = nspace.define(arg)
                    break

            if help_msg is None:
                raise KeyError('invalid argument', arg)

            parser.add_argument(
                '--{}'.format(arg), choices=choices, metavar=get_metavar(value_type),
                help=help_msg, required=required, default=default_value, type=value_type
            )


def generate_config(args, parser, log=devnull):
    """
    Args:
        parser (argparse.ArgumentParser): the main parser
        required: the argparse required arguments group
        optional: the argparse optional arguments group
    """
    libs = []
    inputs_by_lib = {}
    convert = {}
    try:
        # process the libraries by input argument (--input)
        for libconf in [LibraryConfig.parse_args(*a) for a in args.library]:
            if not libconf.bam_file and SUBCOMMAND.VALIDATE not in args.skip_stage:
                raise KeyError('argument --library: bam file must be given if validation is not being skipped')
            libs.append(libconf)
            inputs_by_lib[libconf.library] = set()

        for arg_list in args.input:
            inputfile = arg_list[0]
            for lib in arg_list[1:]:
                if lib not in inputs_by_lib:
                    raise KeyError(
                        'argument --input: specified a library that was not configured. Please input all libraries using '
                        'the --library flag', lib)
                inputs_by_lib[lib].add(inputfile)
        # process the inputs by library argument (--assign)
        for arg_list in args.assign:
            lib = arg_list[0]
            if lib not in inputs_by_lib:
                raise KeyError(
                    'argument --assign: specified a library that was not configured. Please input all libraries using '
                    'the --library flag', lib)
            inputs_by_lib[lib].update(arg_list[1:])

        for libconf in libs:
            if not inputs_by_lib[libconf.library]:
                raise KeyError('argument --input: no input was given for the library', libconf.library)
            libconf.inputs = inputs_by_lib[libconf.library]

        for alias, command in args.external_conversion:
            if alias in convert:
                raise KeyError('duplicate alias names are not allowed', alias)
            convert[alias] = []
            open_option = False
            for item in re.split(r'\s+', command):
                if convert[alias]:
                    if open_option:
                        convert[alias][-1] += ' ' + item
                        open_option = False
                    else:
                        convert[alias].append(item)
                        if item[0] == '-':
                            open_option = True
                else:
                    convert[alias].append(item)

        for arg in args.convert:
            # should follow the pattern: alias file [file...] toolname [stranded]
            alias = arg[0]
            if alias in convert:
                raise KeyError('duplicate alias names are not allowed: {}'.format(alias))
            if arg[-1] in SUPPORTED_TOOL.values():
                toolname = arg[-1]
                stranded = False
                inputfiles = arg[1:-1]
            else:
                toolname, stranded = arg[-2:]
                inputfiles = arg[1:-2]
            if not inputfiles:
                raise KeyError('argument --convert is missing input file path(s): {}'.format(arg))
            stranded = str(tab.cast_boolean(stranded))
            SUPPORTED_TOOL.enforce(toolname)
            convert[alias] = ['convert_tool_output'] + inputfiles + [toolname, stranded]
    except KeyError as err:
        parser.error(' '.join(err.args))

    if SUBCOMMAND.VALIDATE not in args.skip_stage:
        # load the annotations if we need them
        if any([l.is_trans() for l in libs]):
            if not args.get('annotations_filename'):
                parser.error('argument --annotations: is required to gather bam stats for transcriptome libraries')
            log('loading the reference annotations file', args.annotations_filename)
            args.annotations = load_annotations(args.annotations_filename)
        for i, libconf in enumerate(libs):
            log('generating the config section for:', libconf.library)
            libs[i] = LibraryConfig.build(
                library=libconf.library, protocol=libconf.protocol, bam_file=libconf.bam_file,
                inputs=inputs_by_lib[libconf.library], strand_specific=libconf.strand_specific,
                disease_status=libconf.disease_status, annotations=args.annotations, log=log,
                sample_size=args.genome_bins if libconf.protocol == PROTOCOL.GENOME else args.transcriptome_bins,
                distribution_fraction=args.distribution_fraction
            )
    write_config(args.write, include_defaults=args.add_defaults, libraries=libs, conversions=convert, log=log)
