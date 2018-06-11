import argparse
from configparser import ConfigParser, ExtendedInterpolation
from copy import copy as _copy
import logging
import os
import re
import sys
import warnings

import tab

from . import __version__
from .align import SUPPORTED_ALIGNER
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .annotate.file_io import REFERENCE_DEFAULTS
from .bam.cache import BamCache
from .bam.stats import compute_genome_bam_stats, compute_transcriptome_bam_stats
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .constants import DISEASE_STATUS, SUBCOMMAND, PROTOCOL, float_fraction
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .schedule.constants import OPTIONS as SUBMIT_OPTIONS
from .schedule.constants import SCHEDULER
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .tools import SUPPORTED_TOOL
from .util import bash_expands, cast, DEVNULL, MavisNamespace, WeakMavisNamespace, filepath, NullableType
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS


CONVERT_OPTIONS = WeakMavisNamespace()
CONVERT_OPTIONS.add('assume_no_untemplated', True, defn='assume that if not given there is no untemplated sequence between the breakpoints')


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
                    self, 'must have at least {} arguments. Given: {}'.format(self.nmin, values))
        elif not self.nmin <= len(values) <= self.nmax:
            raise argparse.ArgumentError(
                self, 'requires {}-{} arguments. Given: {}'.format(self.nmin, self.nmax, values))
        setattr(namespace, self.dest, items)


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
        self.read_length = NullableType(int)(read_length)
        self.median_fragment_size = NullableType(int)(median_fragment_size)
        self.stdev_fragment_size = NullableType(int)(stdev_fragment_size)
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
                self.add(attr, value, listable=namespace.is_listable(attr), nullable=namespace.is_nullable(attr), cast_type=namespace.type(attr))
                break

    def flatten(self):
        result = MavisNamespace.items(self)
        result['inputs'] = '\n'.join(result['inputs'])
        return result

    def is_trans(self):
        return True if self.protocol == PROTOCOL.TRANS else False

    @staticmethod
    def build(
        library, protocol, bam_file, inputs,
        annotations=None,
        log=DEVNULL,
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

        if protocol == PROTOCOL.TRANS:
            if annotations is None or annotations.is_empty():
                raise AttributeError(
                    'missing required attribute: annotations. Annotations must be given for transcriptomes')
            annotations.load()
        bam = BamCache(bam_file)
        if protocol == PROTOCOL.TRANS:
            bamstats = compute_transcriptome_bam_stats(
                bam,
                annotations=annotations.content,
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


class MavisConfig(MavisNamespace):

    def __init__(self, **kwargs):
        # section can be named schedule or qsub to support older versions
        MavisNamespace.__init__(self)
        try:
            content = validate_section(kwargs.pop('schedule', kwargs.pop('qsub', {})), SUBMIT_OPTIONS, True)
            self.schedule = content
        except Exception as err:
            err.args = ['Error in validating the schedule section in the config. ' + ' '.join([str(a) for a in err.args])]
            raise err

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
            try:
                self[sec] = validate_section(kwargs.pop(sec, {}), defaults, True)
            except Exception as err:
                err.args = ['Error in validating the {} section in the config. '.format(sec) + ' '.join([str(a) for a in err.args])]

                raise err

        SUPPORTED_ALIGNER.enforce(self.validate.aligner)
        for attr, fnames in self.reference.items():
            if attr != 'aligner_reference':
                self.reference[attr] = [f for f in [NullableType(filepath)(v) for v in fnames] if f]
            if not self.reference[attr] and attr not in {'dgv_annotation', 'masking', 'template_metadata'}:
                raise FileNotFoundError(
                    'Error in validating the convert section of the config for tag={}. '
                    'Required reference file does not exist'.format(attr))

        # set the conversion section
        self.convert = kwargs.pop('convert', {})
        for attr, val in self.convert.items():
            if attr in CONVERT_OPTIONS:
                self.convert[attr] = CONVERT_OPTIONS.type(attr)(val)
                continue
            val = [v for v in re.split(r'[;\s]+', val) if v]
            if not val:
                raise UserWarning('Error in validating convert section of the config for tag={}. Tag requires arguments'.format(attr))
            if val[0] == 'convert_tool_output':
                try:
                    val[-1] = tab.cast_boolean(val[-1])
                except TypeError:
                    val.append(False)
                if len(val) < 4 or val[-2] not in SUPPORTED_TOOL.values():
                    raise UserWarning(
                        'Error in validating the convert section of the config for tag={}. '.format(attr),
                        'Conversion using the built-in convert_tool_output requires specifying the input file(s) and '
                        'tool name. Currently supported tools include:', SUPPORTED_TOOL.values(), 'given', val)
                expanded_inputs = []
                for file_expr in val[1:-2]:
                    expanded = bash_expands(file_expr)
                    if not expanded:
                        raise FileNotFoundError(
                            'Error in validating the config for tag={}. '
                            'Input file(s) do not exist'.format(attr), val[1:-2])
                    expanded_inputs.extend(expanded)
                val = [val[0]] + expanded_inputs + val[-2:]
            self.convert[attr] = val
        self.convert = MavisNamespace(**self.convert)

        # now add the library specific sections
        self.libraries = {}

        for libname, val in kwargs.items():  # all other sections already popped
            libname = nameable_string(libname)
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
                    raise UserWarning('Error in validating the library section of the config.', libname, err, terr)
            for inputfile in lc.inputs:
                if inputfile not in self.convert and not os.path.exists(inputfile):
                    raise FileNotFoundError(
                        'Error in validating the library section of the config. Input file does not exist', libname, inputfile)

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
            raise FileNotFoundError('File does not exist: {}'.format(filepath))
        parser = ConfigParser(interpolation=ExtendedInterpolation())
        parser.read(filepath)
        config_dict = {}

        # get the library sections and add the default settings
        for sec in parser.sections():
            config_dict.setdefault(sec, {}).update(parser[sec].items())
        return MavisConfig(**config_dict)


def write_config(filename, include_defaults=False, libraries=[], conversions={}, log=DEVNULL):
    """
    Args:
        filename (str): path to the output file
        include_defaults (bool): True if default parameters should be written to the config, False otherwise
        libraries (list of LibraryConfig): library configuration sections
        conversions (dict of list by str): conversion commands by alias name
        log (function): function to pass output logging to
    """
    config = {}

    config['reference'] = REFERENCE_DEFAULTS.to_dict()
    for filetype, fname in REFERENCE_DEFAULTS.items():
        if fname is None:
            warnings.warn('filetype {} has not been set. This must be done manually before the configuration file is used'.format(filetype))

    if libraries:
        for lib in libraries:
            config[lib.library] = lib.to_dict()

    if include_defaults:
        config['schedule'] = SUBMIT_OPTIONS.to_dict()
        config['validate'] = VALIDATION_DEFAULTS.to_dict()
        config['cluster'] = CLUSTER_DEFAULTS.to_dict()
        config['annotate'] = ANNOTATION_DEFAULTS.to_dict()
        config['illustrate'] = ILLUSTRATION_DEFAULTS.to_dict()
        config['summary'] = SUMMARY_DEFAULTS.to_dict()

    config['convert'] = CONVERT_OPTIONS.to_dict()
    for alias, command in conversions.items():
        if alias in CONVERT_OPTIONS:
            raise UserWarning('error in writing config. Alias for conversion product cannot be a setting', alias, CONVERT_OPTIONS.keys())
        config['convert'][alias] = '\n'.join(command)

    for sec in config:
        for tag, value in config[sec].items():
            if '_regex_' in tag:
                config[sec][tag] = re.sub(r'\$', '$$', config[sec][tag])
                continue
            elif not isinstance(value, str):
                try:
                    config[sec][tag] = '\n'.join([str(v) for v in value])
                    continue
                except TypeError:
                    pass
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
        new_namespace.copy_from(namespace)

    for attr, value in section.items():
        if attr not in namespace:
            raise KeyError('tag not recognized', attr)
        else:
            cast_type = namespace.type(attr)
            if namespace.is_listable(attr):
                value = MavisNamespace.parse_listable_string(value, cast_type, namespace.is_nullable(attr))
            else:
                value = cast_type(value)
            try:
                new_namespace.add(attr, value, cast_type=cast_type, listable=namespace.is_listable(attr), nullable=namespace.is_nullable(attr))
            except Exception as err:
                raise ValueError('failed adding {}. {}'.format(attr, err))
    return new_namespace


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


def nameable_string(input_string):
    """
    A string that can be used for library and/or filenames
    """
    input_string = str(input_string)
    if re.search(r'[;,_\s]', input_string):
        raise TypeError('names cannot contain the reserved characters [;,_\\s]', input_string)
    if input_string.lower() == 'none':
        raise TypeError('names cannot be none', input_string)
    if not input_string:
        raise TypeError('names cannot be an empty string', input_string)
    if not re.search(r'^[a-zA-Z]', input_string):
        raise TypeError('names must start with a letter', input_string)
    return input_string


def augment_parser(arguments, parser, required=None):
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
        elif arg == 'log':
            parser.add_argument('--log', help='redirect stdout to a log file', default=None)
        elif arg == 'log_level':
            parser.add_argument('--log_level', help='level of logging to output', choices=['INFO', 'DEBUG'], default='INFO')
        elif arg == 'aligner_reference':
            default = REFERENCE_DEFAULTS[arg]
            parser.add_argument(
                '--{}'.format(arg), default=default, required=required if not default else False,
                help=REFERENCE_DEFAULTS.define(arg), type=filepath)
        elif arg in REFERENCE_DEFAULTS:
            default = REFERENCE_DEFAULTS[arg]
            parser.add_argument(
                '--{}'.format(arg), default=default, required=required if not default else False,
                help=REFERENCE_DEFAULTS.define(arg), type=filepath if required else NullableType(filepath), nargs='*')
        elif arg == 'config':
            parser.add_argument('config', help='path to the config file', type=filepath)
        elif arg == 'bam_file':
            parser.add_argument('--bam_file', help='path to the input bam file', required=required, type=filepath)
        elif arg == 'read_length':
            parser.add_argument(
                '--read_length', type=int, help='the length of the reads in the bam file',
                required=required)
        elif arg == 'stdev_fragment_size':
            parser.add_argument(
                '--stdev_fragment_size', type=int, help='expected standard deviation in insert sizes',
                required=required)
        elif arg == 'median_fragment_size':
            parser.add_argument(
                '--median_fragment_size', type=int, help='median inset size for pairs in the bam file', required=required)
        elif arg == 'library':
            parser.add_argument('--library', help='library name', required=required, type=nameable_string)
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
                '--strand_specific', type=tab.cast_boolean,
                default=False, help='indicates that the input is strand specific')
        else:
            value_type = None
            help_msg = None
            default_value = None
            choices = None
            nargs = None
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
                    if nspace.is_listable(arg):
                        nargs = '*'
                    value_type = nspace.type(arg, None)
                    if nspace.is_nullable(arg):
                        value_type = NullableType(value_type)
                    if not help_msg:
                        help_msg = nspace.define(arg)
                    break

            if help_msg is None:
                raise KeyError('invalid argument', arg)
            parser.add_argument(
                '--{}'.format(arg), choices=choices, nargs=nargs,
                help=help_msg, required=required, default=default_value, type=value_type
            )


def generate_config(args, parser, log=DEVNULL):
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
            if SUBCOMMAND.VALIDATE not in args.skip_stage and libconf.protocol == PROTOCOL.TRANS and (not args.annotations or args.annotations.is_empty()):
                parser.error('argument --annotations is required to build configuration files for transcriptome libraries')

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
