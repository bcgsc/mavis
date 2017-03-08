from argparse import Namespace
from configparser import ConfigParser, ExtendedInterpolation
import argparse
import os
import TSV
from .. import __version__
from ..constants import PROTOCOL
from ..validate.constants import VALIDATION_DEFAULTS
from .util import get_blat_version, get_samtools_version, PIPELINE_STEP

QSUB_TAGS = dict(memory=12, queue='transabyss.q')

LIBRARY_DEFAULT_TAGS = dict(
    min_clusters_per_file=50,
    max_files=10,
    cluster_clique_size=15,
    cluster_radius=20,
    min_orf_size=120,
    max_orf_cap=3,
    min_domain_mapping_match=0.8,
    max_proximity=5000,
    uninformative_filter=True,
    stranded_bam=False,
    domain_regex_filter='^PF\d+$'
)
LIBRARY_DEFAULT_TAGS.update(VALIDATION_DEFAULTS.__dict__)

REFERENCE_TAGS = ['template_metadata', 'reference_genome', 'annotations', 'masking', 'blat_2bit_reference']

LIBRARY_REQUIRED_TAGS = dict(
    protocol=PROTOCOL.enforce,
    bam_file=str, 
    read_length=int, 
    median_fragment_size=int, 
    stdev_fragment_size=int, 
    inputs=lambda x: x.split(';') if x else [], 
    pairing=lambda x: x.split(';') if x else []
)


def write_config(filename, include_defaults=False):
    config = ConfigParser()
    
    for sec in ['DEFAULTS', 'reference', '<LIBRARY NAME>', 'qsub']:
        config[sec] = {}
    
    for tag in REFERENCE_TAGS:
        config['reference'][tag] = '<REQUIRED>'
    for tag in LIBRARY_REQUIRED_TAGS:
        config['<LIBRARY NAME>'][tag] = '<REQUIRED>'
    
    if include_defaults:
        for tag, val in QSUB_TAGS.items():
            config['qsub'][tag] = str(val)
        for tag, val in LIBRARY_DEFAULT_TAGS.items():
            config['DEFAULTS'][tag] = str(val)
    
    with open(filename, 'w') as configfile:
        config.write(configfile)


def cast(value, cast_func):
    if cast_func == bool:
        value = TSV.tsv_boolean(value)
    else:
        value = cast_func(value)
    return value


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
        if sec not in ['DEFAULTS', 'reference', 'qsub']:
            library_sections.append(sec)
    
    all_libs = {}
    all_libs.update(QSUB_TAGS)
    all_libs.update(LIBRARY_DEFAULT_TAGS)
    # check that the reference files all exist
    for attr, fname in parser['reference'].items():
        if not os.path.exists(fname):
            raise KeyError(attr, 'file at', fname, 'dose not exist')
        all_libs[attr] = fname
    for attr in REFERENCE_TAGS:
        if attr not in parser['reference']:
            raise KeyError('missing required tag', attr, 'in reference section')

    # type check the qsub options
    for attr, value in parser['qsub'].items():
        if attr not in QSUB_TAGS:
            raise KeyError('unrecognized tag in the sub section', attr)
        all_libs[attr] = cast(value, type(QSUB_TAGS[attr]))
    
    # cast the defaults
    for attr, value in parser['DEFAULTS'].items():
        if attr not in LIBRARY_DEFAULT_TAGS:
            raise KeyError('unrecognized tag in the DEFAULTS section', attr)
        all_libs[attr] = cast(value, type(LIBRARY_DEFAULT_TAGS[attr]))
    
    sections = []
    for sec in library_sections:
        d = {}
        d.update(all_libs)
        for attr, value in parser[sec].items():
            if attr in LIBRARY_DEFAULT_TAGS:
                value = cast(value, type(LIBRARY_DEFAULT_TAGS[attr]))
            d[attr] = value
        for attr in LIBRARY_REQUIRED_TAGS:
            if attr not in parser[sec]:
                raise KeyError('required tag', attr, 'not found in library section', sec)
            else:
                d[attr] = LIBRARY_REQUIRED_TAGS[attr](parser[sec][attr])
        d['library'] = sec
        sections.append(Namespace(**d))
    if len(library_sections) < 1:
        raise UserWarning('configuration file must have 1 or more library sections')
    return sections


def parse_arguments(pstep):
    PIPELINE_STEP.enforce(pstep)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )

    if pstep != PIPELINE_STEP.PIPELINE:
        g = parser.add_argument_group('reference input arguments')
        g.add_argument(
            '--annotations',
            default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
            help='path to the reference annotations of genes, transcript, exons, domains, etc.'
        )
        if pstep in [PIPELINE_STEP.ANNOTATE, PIPELINE_STEP.VALIDATE]:
            g.add_argument(
                '--reference_genome',
                default='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
                help='path to the human reference genome in fa format'
            )
        if pstep == PIPELINE_STEP.ANNOTATE:
            g.add_argument(
                '--template_metadata', default=os.path.join(os.path.dirname(__file__), 'cytoBand.txt'),
                help='file containing the cytoband template information'
            )
        if pstep in [PIPELINE_STEP.CLUSTER, PIPELINE_STEP.VALIDATE]:
            g.add_argument(
                '--masking',
                default='/home/creisle/svn/svmerge/trunk/hg19_masked_regions.tsv'
            )
        g.add_argument(
            '--low_memory', default=False, type=TSV.tsv_boolean,
            help='when working on a machine with less memory this is sacrifice time for memory where possible'
        )
        if pstep == PIPELINE_STEP.VALIDATE:
            g.add_argument(
                '--blat_2bit_reference', default='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
                help='path to the 2bit reference file used for blatting contig sequences'
            )
    else:
        parser.add_argument('config', help='path to the pipeline configuration file')
        parser.add_argument(
            '-f', '--force_overwrite', default=False, type=TSV.tsv_boolean,
            help='set flag to overwrite existing reviewed files'
        )
    if pstep == PIPELINE_STEP.PIPELINE:
        m = parser.add_mutually_exclusive_group(required=True)
        m.add_argument('--output', help='path to the output directory')
        m.add_argument('--write', default=False, action='store_true', help='write a config')
    else:
        parser.add_argument('--output', help='path to the output directory', required=True)

    if pstep == PIPELINE_STEP.ANNOTATE:
        parser.add_argument(
            '--output_svgs', default=True, type=TSV.tsv_boolean,
            help='set flag to suppress svg drawings of putative annotations')
        parser.add_argument(
            '--min_orf_size', default=LIBRARY_DEFAULT_TAGS['min_orf_size'], type=int, 
            help='minimum size for putative ORFs'
        )
        parser.add_argument(
            '--max_orf_cap', default=LIBRARY_DEFAULT_TAGS['max_orf_cap'], type=int, 
            help='keep the n longest orfs'
        )
        parser.add_argument(
            '--min_domain_mapping_match', default=LIBRARY_DEFAULT_TAGS['min_domain_mapping_match'], type=float,
            help='minimum percent match for the domain to be considered aligned'
        )
        g = parser.add_argument_group('visualization options')
        g.add_argument(
            '--domain_regex_filter', default='^PF\d+$',
            help='only show domains which names (external identifiers) match the given pattern'
        )

    if pstep == PIPELINE_STEP.CLUSTER or pstep == PIPELINE_STEP.ANNOTATE or pstep == PIPELINE_STEP.PAIR:
        parser.add_argument(
            '--max_proximity', default=LIBRARY_DEFAULT_TAGS['max_proximity'], type=int,
            help='The maximum distance away from breakpoints to look for proximal genes'
        )
        parser.add_argument('-n', '--inputs', required=True, nargs='+', help='1 or more input files')
    elif pstep == PIPELINE_STEP.VALIDATE:
        parser.add_argument('-n', '--input', help='path to the input file', required=True)

    if pstep == PIPELINE_STEP.CLUSTER or pstep == PIPELINE_STEP.VALIDATE:
        parser.add_argument('-l', '--library', help='library name')
        parser.add_argument('--protocol', help='the library protocol: genome or transcriptome', choices=PROTOCOL.values())

    if pstep == PIPELINE_STEP.CLUSTER:
        g = parser.add_argument_group('output arguments')
        g.add_argument(
            '--max_files', default=LIBRARY_DEFAULT_TAGS['max_files'], type=int, dest='max_files',
            help='defines the maximum number of files that can be created')
        g.add_argument(
            '--min_clusters_per_file', default=LIBRARY_DEFAULT_TAGS['min_clusters_per_file'], type=int,
            help='defines the minimum number of clusters per file')
        parser.add_argument(
            '-r', '--cluster_radius', help='radius to use in clustering', 
            default=LIBRARY_DEFAULT_TAGS['cluster_radius'], type=int)
        parser.add_argument(
            '-k', '--cluster_clique_size',
            help='parameter used for computing cliques, smaller is faster, above 20 will be slow',
            default=LIBRARY_DEFAULT_TAGS['cluster_clique_size'], type=int
        )
        parser.add_argument(
            '--uninformative_filter', default=LIBRARY_DEFAULT_TAGS['uninformative_filter'], 
            help='If flag is False then the clusters will not be filtered '
            'based on lack of annotation', type=TSV.tsv_boolean
        )

    if pstep == PIPELINE_STEP.PAIR:
        parser.add_argument(
            '--split_call_distance', default=10, type=int,
            help='distance allowed between breakpoint calls when pairing from split read (and higher) resolution calls'
        )
        parser.add_argument(
            '--contig_call_distance', default=0, type=int,
            help='distance allowed between breakpoint calls when pairing from contig (and higher) resolution calls'
        )
        parser.add_argument(
            '--flanking_call_distance', default=0, type=int,
            help='distance allowed between breakpoint calls when pairing from contig (and higher) resolution calls'
        )

    if pstep == PIPELINE_STEP.VALIDATE:
        g = parser.add_argument_group('evidence arguments')
        for attr, value in VALIDATION_DEFAULTS.__dict__.items():
            vtype = type(value)
            if type(value) == bool:
                vtype = TSV.tsv_boolean
            g.add_argument('--{}'.format(attr), default=value, type=vtype, help='see user manual for desc')
        parser.add_argument(
            '-b', '--bam_file',
            help='path to the input bam file', required=True
        )
        parser.add_argument(
            '--stranded_bam', default=False, type=TSV.tsv_boolean,
            help='indicates that the input bam file is strand specific'
        )
        g.add_argument('--read_length', type=int, help='the length of the reads in the bam file', required=True)
        g.add_argument(
            '--stdev_fragment_size', type=int, help='expected standard deviation in insert sizes', required=True
        )
        g.add_argument(
            '--median_fragment_size', type=int, help='median inset size for pairs in the bam file', required=True
        )

    args = parser.parse_args()
    if pstep == PIPELINE_STEP.VALIDATE:
        args.samtools_version = get_samtools_version()
        args.blat_version = get_blat_version()
    try:
        args.output = os.path.abspath(args.output)
        if os.path.exists(args.output) and not args.force_overwrite:
            parser.print_help()
            print(
                '\nerror: output directory {} exists, --force_overwrite must be specified or the directory removed'.format(
                    repr(args.output)))
            exit(1)
    except (AttributeError, TypeError):
        pass

    return args


