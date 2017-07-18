from configparser import ConfigParser, ExtendedInterpolation
import os
import TSV
import re
import pysam
from . import __version__
from .constants import PROTOCOL, DISEASE_STATUS
from .util import devnull, MavisNamespace
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .tools import SUPPORTED_TOOL
from .align import SUPPORTED_ALIGNER

from .bam.stats import compute_genome_bam_stats, compute_transcriptome_bam_stats

ENV_VAR_PREFIX = 'MAVIS_'


def cast(value, cast_func):
    if cast_func == bool:
        value = TSV.tsv_boolean(value)
    else:
        value = cast_func(value)
    return value


class LibraryConfig:
    def __init__(
        self, library, protocol, disease_status, bam_file, inputs, read_length, median_fragment_size,
        stdev_fragment_size, stranded_bam, strand_determining_read=2,
        **kwargs
    ):
        self.library = library
        self.protocol = PROTOCOL.enforce(protocol)
        self.bam_file = bam_file
        self.read_length = int(read_length)
        self.median_fragment_size = int(median_fragment_size)
        self.stdev_fragment_size = int(stdev_fragment_size)
        self.stranded_bam = cast(stranded_bam, bool)
        self.strand_determining_read = int(strand_determining_read)
        self.disease_status = DISEASE_STATUS.enforce(disease_status)
        try:
            self.inputs = [f for f in re.split('[;\s]+', inputs) if f]
        except TypeError:
            self.inputs = inputs

        acceptable = {}
        acceptable.update(VALIDATION_DEFAULTS.__dict__)

        for attr, value in kwargs.items():
            if attr == 'assembly_max_kmer_size' and value in [None, 'None', '']:  # special case
                setattr(self, attr, None)
            else:
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

    @classmethod
    def build(
        self, library, protocol, bam_file, inputs,
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
        try:
            bam = pysam.AlignmentFile(bam_file, 'rb')
            bamstats = None
            if protocol == PROTOCOL.TRANS:
                bamstats = compute_transcriptome_bam_stats(
                    bam,
                    annotations=annotations,
                    sample_size=sample_size,
                    sample_cap=sample_cap,
                    distribution_fraction=distribution_fraction,
                    log=log
                )
            elif protocol == PROTOCOL.GENOME:
                bamstats = compute_genome_bam_stats(
                    bam,
                    sample_size=sample_size,
                    sample_bin_size=sample_bin_size,
                    sample_cap=sample_cap,
                    distribution_fraction=distribution_fraction,
                    log=log
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
        finally:
            try:
                bam.close()
            except AttributeError:
                pass


class JobSchedulingConfig:
    def __init__(self, validate_memory_gb=12, default_memory_gb=10, queue='transabyss.q'):
        self.validate_memory_gb = validate_memory_gb
        self.default_memory_gb = default_memory_gb
        self.queue = queue

    def flatten(self):
        result = {}
        result.update(self.__dict__)
        return result


class ReferenceFilesConfig:

    def __init__(
        self,
        annotations=None,
        reference_genome=None,
        template_metadata=None,
        masking=None,
        aligner_reference=None,
        dgv_annotation=None,
        low_memory=False
    ):
        self.annotations = annotations or os.environ.get(ENV_VAR_PREFIX + 'ANNOTATIONS', None)
        self.reference_genome = reference_genome or os.environ.get(ENV_VAR_PREFIX + 'REFERENCE_GENOME', None)
        self.template_metadata = template_metadata or os.environ.get(ENV_VAR_PREFIX + 'TEMPLATE_METADATA', None)
        self.masking = masking or os.environ.get(ENV_VAR_PREFIX + 'MASKING', None)
        self.dgv_annotation = dgv_annotation or os.environ.get(ENV_VAR_PREFIX + 'DGV_ANNOTATION', None)
        self.low_memory = low_memory or os.environ.get(ENV_VAR_PREFIX + 'LOW_MEMORY', None)
        self.aligner_reference = aligner_reference or os.environ.get(ENV_VAR_PREFIX + 'ALIGNER_REFERENCE', None)

    def flatten(self):
        result = {}
        result.update(self.__dict__)
        return result


class PairingConfig:

    def __init__(
        self,
        split_call_distance=PAIRING_DEFAULTS.split_call_distance,
        contig_call_distance=PAIRING_DEFAULTS.contig_call_distance,
        flanking_call_distance=PAIRING_DEFAULTS.flanking_call_distance,
        spanning_call_distance=PAIRING_DEFAULTS.spanning_call_distance,
        max_proximity=CLUSTER_DEFAULTS.max_proximity,
        low_memory=False
    ):
        self.split_call_distance = int(split_call_distance)
        self.contig_call_distance = int(contig_call_distance)
        self.flanking_call_distance = int(flanking_call_distance)
        self.spanning_call_distance = int(spanning_call_distance)

    def flatten(self):
        result = {}
        result.update(self.__dict__)
        return result


class SummaryConfig:
    def __init__(
        self,
        filter_min_remapped_reads=SUMMARY_DEFAULTS.filter_min_remapped_reads,
        filter_min_spanning_reads=SUMMARY_DEFAULTS.filter_min_spanning_reads,
        filter_min_flanking_reads=SUMMARY_DEFAULTS.filter_min_flanking_reads,
        filter_min_flanking_only_reads=SUMMARY_DEFAULTS.filter_min_flanking_only_reads,
        filter_min_split_reads=SUMMARY_DEFAULTS.filter_min_split_reads,
        filter_min_linking_split_reads=SUMMARY_DEFAULTS.filter_min_linking_split_reads,
        flanking_call_distance=PAIRING_DEFAULTS.flanking_call_distance,
        split_call_distance=PAIRING_DEFAULTS.split_call_distance,
        contig_call_distance=PAIRING_DEFAULTS.contig_call_distance,
        spanning_call_distance=PAIRING_DEFAULTS.spanning_call_distance,
    ):
        self.filter_min_remapped_reads = int(filter_min_remapped_reads)
        self.filter_min_spanning_reads = int(filter_min_spanning_reads)
        self.filter_min_flanking_reads = int(filter_min_flanking_reads)
        self.filter_min_flanking_only_reads = int(filter_min_flanking_only_reads)
        self.filter_min_split_reads = int(filter_min_split_reads)
        self.filter_min_linking_split_reads = int(filter_min_linking_split_reads)

    def flatten(self):
        result = {}
        result.update(self.__dict__)
        return result


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

    config['reference'] = ReferenceFilesConfig().flatten()

    if libraries:
        for lib in libraries:
            config[lib.library] = lib.flatten()

    if include_defaults:
        config['qsub'] = JobSchedulingConfig().flatten()
        config['validation'] = {}
        config['validation'].update(VALIDATION_DEFAULTS.__dict__)
        config['cluster'] = {}
        config['cluster'].update(CLUSTER_DEFAULTS.__dict__)
        config['annotation'] = {}
        config['annotation'].update(ANNOTATION_DEFAULTS.__dict__)
        config['illustrate'] = {}
        config['illustrate'].update(ILLUSTRATION_DEFAULTS.__dict__)
        for sec in ['qsub', 'illustrate', 'validation', 'cluster', 'annotation']:
            for tag, val in config[sec].items():
                env = ENV_VAR_PREFIX + tag.upper()
                config[sec][tag] = os.environ.get(env, None) or val
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

    with open(filename, 'w') as configfile:
        log('writing:', filename)
        conf.write(configfile)


def validate_and_cast_section(section, defaults):
    d = {}
    for attr, value in section.items():
        if attr not in defaults:
            raise KeyError('tag not recognized', attr)
        else:
            d[attr] = cast(value, type(defaults[attr]))
    return d


def read_config(filepath):
    """
    reads the configuration settings from the configuration file

    Args:
        filepath (str): path to the input configuration file

    Returns:
        class:`list` of :class:`Namespace`: namespace arguments for each library
    """
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read(filepath)

    # get the library sections and add the default settings
    library_sections = []
    for sec in parser.sections():
        if sec not in [
            'validation', 'reference', 'qsub', 'illustrate',
            'annotation', 'cluster', 'annotation', 'convert'
        ]:
            library_sections.append(sec)
            parser[sec]['library'] = sec

    job_sched = JobSchedulingConfig(**(parser['qsub'] if 'qsub' in parser else {}))
    ref = ReferenceFilesConfig(**(parser['reference'] if 'reference' in parser else {}))
    pairing = PairingConfig(**(parser['pairing'] if 'pairing' in parser else {}))
    summary = SummaryConfig(**(parser['summary'] if 'summary' in parser else {}))
    global_args = {}
    global_args.update(job_sched.flatten())
    global_args.update(ref.flatten())
    global_args.update(ILLUSTRATION_DEFAULTS.__dict__)
    global_args.update(pairing.flatten())
    global_args.update(summary.flatten())
    global_args.update(ANNOTATION_DEFAULTS.__dict__)
    global_args.update(CLUSTER_DEFAULTS.__dict__)
    try:
        global_args.update(validate_and_cast_section(parser['illustrate'], ILLUSTRATION_DEFAULTS))
    except KeyError:
        pass

    try:
        global_args.update(validate_and_cast_section(parser['annotation'], ANNOTATION_DEFAULTS))
    except KeyError:
        pass

    try:
        global_args.update(validate_and_cast_section(parser['cluster'], CLUSTER_DEFAULTS))
    except KeyError:
        pass

    args = {}
    args.update(VALIDATION_DEFAULTS.__dict__)
    try:
        args.update(parser['validation'] if 'validation' in parser else {})
        SUPPORTED_ALIGNER.enforce(args['aligner'])
    except KeyError:
        pass

    # check that the reference files all exist
    for attr, fname in parser['reference'].items():
        if not os.path.exists(fname) and attr != 'low_memory':
            raise KeyError(attr, 'file at', fname, 'does not exist')
        global_args[attr] = fname

    convert = {}
    if 'convert' in parser:
        for attr, val in parser['convert'].items():
            val = [v for v in re.split('[;\s]+', val) if v]
            if val[0] == 'convert_tool_output':
                if len(val) < 3 or val[2] not in SUPPORTED_TOOL:
                    raise UserWarning(
                        'conversion using the built-in convert_tool_output requires specifying the input file and '
                        'tool name currently supported tools include:', SUPPORTED_TOOL.values())
                elif len(val) == 4:
                    val[3] = TSV.tsv_boolean(val[3])
                else:
                    raise UserWarning(
                        'conversion using the built-in convert_tool_output takes at most 3 arguments')
            convert[attr] = val

    sections = []
    for sec in library_sections:
        d = {}
        d.update(args)
        d.update(parser[sec])

        # now try building the LibraryConfig object
        try:
            lc = LibraryConfig(**d)
            sections.append(lc)
            continue
        except TypeError as terr:  # missing required argument
            try:
                lc = LibraryConfig.build(**d)
                sections.append(lc)
            except Exception as err:
                raise UserWarning('could not build configuration file', terr, err)
    for sec in sections:
        for infile in sec.inputs:
            if not os.path.exists(infile) and infile not in convert:
                raise UserWarning('input file does not exist and is not a conversion', infile)
    if len(library_sections) < 1:
        raise UserWarning('configuration file must have 1 or more library sections')
    return MavisNamespace(**global_args), sections, MavisNamespace(**convert)


def add_semi_optional_argument(argname, success_parser, failure_parser, help_msg=''):
    """
    for an argument tries to get the argument default from the environment variable
    """
    env_name = ENV_VAR_PREFIX + argname.upper()
    help_msg += ' The default for this argument is configured by setting the environment variable {}'.format(env_name)
    if os.environ.get(env_name, None):
        success_parser.add_argument('--{}'.format(argname), required=False, default=os.environ[env_name], help=help_msg)
    else:
        failure_parser.add_argument('--{}'.format(argname), required=True, help=help_msg)


def get_env_variable(arg, default, cast_type=None):
    """
    Args:
        arg (str): the argument/variable name
    Returns:
        the setting from the environment variable if given, otherwise the default value
    """
    if cast_type is None:
        cast_type = type(default)
    name = ENV_VAR_PREFIX + arg.upper()
    result = os.environ.get(name, None)
    if result is not None:
        return cast(result, cast_type)
    else:
        return default


def augment_parser(parser, optparser, arguments):
    for arg in arguments:
        if arg == 'help':
            optparser.add_argument('-h', '--help', action='help', help='show this help message and exit')
        elif arg == 'version':
            optparser.add_argument(
                '-v', '--version', action='version', version='%(prog)s version ' + __version__,
                help='Outputs the version number')
        elif arg == 'annotations':
            add_semi_optional_argument(
                arg, optparser, parser, 'Path to the reference annotations of genes, transcript, exons, domains, etc.')
        elif arg == 'reference_genome':
            add_semi_optional_argument(arg, optparser, parser, 'Path to the human reference genome fasta file.')
            optparser.add_argument(
                '--low_memory', default=get_env_variable('low_memory', False), type=TSV.tsv_boolean,
                help='if true defaults to indexing vs loading the reference genome')
        elif arg == 'template_metadata':
            add_semi_optional_argument(arg, optparser, parser, 'File containing the cytoband template information.')
        elif arg == 'masking':
            add_semi_optional_argument(arg, optparser, parser)
        elif arg == 'aligner':
            optparser.add_argument(
                '--' + arg, default=get_env_variable(arg, VALIDATION_DEFAULTS[arg]),
                choices=SUPPORTED_ALIGNER.values(), help='aligner to use for aligning contigs')
        elif arg == 'aligner_reference':
            add_semi_optional_argument(
                arg, optparser, parser, 'path to the aligner reference file used for aligning the contig sequences.')
        elif arg == 'dgv_annotation':
            add_semi_optional_argument(
                arg, optparser, parser, 'Path to the dgv reference processed to look like the cytoband file.')
        elif arg == 'config':
            parser.add_argument('config', 'path to the config file')
        elif arg == 'stranded_bam':
            optparser.add_argument(
                '--stranded_bam', required=True, type=TSV.tsv_boolean,
                help='indicates that the input bam file is strand specific')
        elif arg == 'force_overwrite':
            optparser.add_argument(
                '-f', '--force_overwrite', default=get_env_variable(arg, False), type=TSV.tsv_boolean,
                help='set flag to overwrite existing reviewed files')
        elif arg == 'output_svgs':
            optparser.add_argument(
                '--output_svgs', default=get_env_variable(arg, True), type=TSV.tsv_boolean,
                help='set flag to suppress svg drawings of putative annotations')
        elif arg == 'max_files':
            optparser.add_argument(
                '--max_files', default=get_env_variable(arg, CLUSTER_DEFAULTS.max_files), type=int, dest='max_files',
                help='defines the maximum number of files that can be created')
        elif arg == 'min_clusters_per_file':
            optparser.add_argument(
                '--min_clusters_per_file', default=get_env_variable(arg, CLUSTER_DEFAULTS.min_clusters_per_file),
                type=int, help='defines the minimum number of clusters per file')
        elif arg == 'cluster_radius':
            optparser.add_argument(
                '-r', '--cluster_radius', help='radius to use in clustering',
                default=get_env_variable(arg, CLUSTER_DEFAULTS.cluster_radius), type=int)
        elif arg == 'cluster_clique_size':
            optparser.add_argument(
                '-k', '--cluster_clique_size', default=get_env_variable(arg, CLUSTER_DEFAULTS.cluster_clique_size),
                type=int, help='parameter used for computing cliques, smaller is faster, above 20 will be slow')
        elif arg == 'uninformative_filter':
            optparser.add_argument(
                '--uninformative_filter', default=get_env_variable(arg, CLUSTER_DEFAULTS.uninformative_filter),
                type=TSV.tsv_boolean,
                help='If flag is False then the clusters will not be filtered based on lack of annotation'
            )
        elif arg in PAIRING_DEFAULTS:
            optparser.add_argument(
                '--{}'.format(arg), default=get_env_variable(arg, PAIRING_DEFAULTS[arg]), type=int,
                help='distance allowed between breakpoint calls when pairing breakpoints of this call method')
        elif arg in VALIDATION_DEFAULTS:
            value = VALIDATION_DEFAULTS[arg]
            vtype = type(value) if not isinstance(value, bool) else TSV.tsv_boolean
            optparser.add_argument(
                '--{}'.format(arg), default=get_env_variable(arg, value), type=vtype, help='see user manual for desc')
        elif arg in ILLUSTRATION_DEFAULTS:
            value = ILLUSTRATION_DEFAULTS[arg]
            vtype = type(value) if not isinstance(value, bool) else TSV.tsv_boolean
            optparser.add_argument(
                '--{}'.format(arg), default=get_env_variable(arg, value), type=vtype, help='see user manual for desc')
        elif arg in ANNOTATION_DEFAULTS:
            value = ANNOTATION_DEFAULTS[arg]
            vtype = type(value) if not isinstance(value, bool) else TSV.tsv_boolean
            optparser.add_argument(
                '--{}'.format(arg), default=get_env_variable(arg, value), type=vtype, help='see user manual for desc')
        elif arg in SUMMARY_DEFAULTS:
            value = SUMMARY_DEFAULTS[arg]
            vtype = type(value) if not isinstance(value, bool) else TSV.tsv_boolean
            optparser.add_argument(
                '--{}'.format(arg), default=get_env_variable(arg, value), type=vtype, help='see user manual for desc')
        elif arg == 'max_proximity':
            optparser.add_argument(
                '--{}'.format(arg), default=get_env_variable(arg, CLUSTER_DEFAULTS[arg]), type=int,
                help='maximum distance away from an annotation before the uninformative filter is applied or the'
                'annotation is not considered for a given event')
        elif arg == 'bam_file':
            parser.add_argument('--bam_file', help='path to the input bam file', required=True)
        elif arg == 'read_length':
            parser.add_argument(
                '--read_length', type=int, help='the length of the reads in the bam file', required=True)
        elif arg == 'stdev_fragment_size':
            parser.add_argument(
                '--stdev_fragment_size', type=int, help='expected standard deviation in insert sizes', required=True)
        elif arg == 'median_fragment_size':
            parser.add_argument(
                '--median_fragment_size', type=int, help='median inset size for pairs in the bam file', required=True)
        elif arg == 'library':
            parser.add_argument('--library', help='library name', required=True)
        elif arg == 'protocol':
            parser.add_argument('--protocol', choices=PROTOCOL.values(), help='library protocol', required=True)
        elif arg == 'disease_status':
            parser.add_argument(
                '--disease_status', choices=DISEASE_STATUS.values(), help='library disease status', required=True)
        else:
            raise KeyError('invalid argument', arg)
