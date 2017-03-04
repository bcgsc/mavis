import argparse
import os
import sys
from configparser import ConfigParser, ExtendedInterpolation
import re
import warnings
from argparse import Namespace
import subprocess
from datetime import datetime

import TSV
from vocab import Vocab

# local modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis.breakpoint import read_bpp_from_input_file, BreakpointPair
from mavis.annotate import load_reference_genes, load_reference_genome, load_masking_regions, load_templates
from mavis import __version__
from mavis.constants import PROTOCOL, COLUMNS, sort_columns
from mavis.validate.constants import VALIDATION_DEFAULTS
from mavis.interval import Interval



MIN_CLUSTERS_PER_FILE = 50
MAX_FILES = 100
CLUSTER_CLIQUE_SIZE = 15
CLUSTER_RADIUS = 20
MIN_ORF_SIZE = 120
MAX_ORF_CAP = 3
MAX_PROXIMITY = 5000
DEFAULTS = Namespace(
    min_clusters_per_file=50,
    max_files=10,
    cluster_clique_size=15,
    cluster_radius=20,
    min_orf_size=120,
    max_orf_cap=3,
    min_domain_mapping_match=0.8,
    domain_regex_filter='^PF\d+$',
    max_proximity=5000,
    uninformative_filter=True,
    blat_2bit_reference='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
    stranded_bam=False
)
DEFAULTS.__dict__.update(VALIDATION_DEFAULTS.__dict__)
OUTPUT_FILES = Namespace(
    CLUSTER_SET='{prefix}'
)

PIPELINE_STEP = Vocab(
    ANNOTATE='annotate',
    VALIDATE='validate',
    PIPELINE='pipeline',
    CLUSTER='cluster',
    PAIR='pairing',
    SUMMARY='summary'
)


def build_batch_id(prefix='', suffix='', size=6):
    date = datetime.now()
    m = int(math.pow(10, size) - 1)
    return 'batch{prefix}{date.year}{date.month:02d}{date.day:02d}r{r:06d}{suffix}'.format(
        prefix=prefix, suffix=suffix, date=date, r=random.randint(1, m))


def samtools_v0_sort(input_bam, output_bam):
    prefix = re.sub('\.bam$', '', output_bam)
    return 'samtools sort {} {}'.format(input_bam, prefix)


def samtools_v1_sort(input_bam, output_bam):
    return 'samtools sort {} -o {}'.format(input_bam, output_bam)


def get_samtools_version():
    proc = subprocess.getoutput(['samtools'])
    for line in proc.split('\n'):
        m = re.search('Version: (?P<major>\d+)\.(?P<mid>\d+)\.(?P<minor>\d+)', line)
        if m:
            return int(m.group('major')), int(m.group('mid')), int(m.group('minor'))
    raise ValueError('unable to parse samtools version number')


def get_blat_version():
    proc = subprocess.getoutput(['blat'])
    for line in proc.split('\n'):
        m = re.search('blat - Standalone BLAT v. (\d+(x\d+)?)', line)
        if m:
            return m.group(1)
    raise ValueError('unable to parse blat version number')


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def mkdirp(dirname):
    log("creating output directory: '{}'".format(dirname))
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise exc


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

    LIBRARY_REQ_ATTR = ['protocol', 'bam_file', 'read_length', 'median_fragment_size', 'stdev_fragment_size', 'inputs']
    TYPE_CHECK = DEFAULTS.__dict__

    config = {
        'qsub': {
            'memory': 12,
            'queue': 'transabyss.q'
        },
        'reference': {
            'template_metadata': os.path.join(os.path.dirname(__file__), 'cytoBand.txt'),
            'reference_genome': '/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
            'annotations': '/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
            'masking': '/home/creisle/svn/svmerge/trunk/hg19_masked_regions.tsv'
        }
    }

    defaults = dict()
    defaults.update(DEFAULTS.__dict__)
    for attr, value in parser['DEFAULTS'].items():
        if attr == 'protocol':
            PROTOCOL.enforce(value)
        if attr in TYPE_CHECK and type(TYPE_CHECK[attr]) != type(value):
            try:
                if type(TYPE_CHECK[attr]) == bool:
                    value = TSV.tsv_boolean(value)
                else:
                    value = type(TYPE_CHECK[attr])(value)
            except ValueError:
                warnings.warn('type check failed for attr {} with value {}'.format(attr, repr(value)))
        elif attr not in TYPE_CHECK:
            raise ValueError('unexpected value in DEFAULTS section', attr, value)
        defaults[attr] = value

    library_sections = []

    for sec in parser.sections():

        section = dict()
        if sec == 'DEFAULTS':
            continue
        elif sec not in ['reference', 'qsub', 'visualization']:  # assume this is a library configuration
            library_sections.append(sec)
            for attr in LIBRARY_REQ_ATTR:
                if not parser.has_option(sec, attr):
                    raise KeyError(
                        'missing one or more required attribute(s) for the library section',
                        sec, attr, LIBRARY_REQ_ATTR)
                if attr in config['reference']:
                    raise ValueError(
                        'this attribute cannot be given per library and must be general between libraries', attr)
            section['library'] = sec

        for attr, value in parser[sec].items():
            if attr == 'protocol':
                PROTOCOL.enforce(value)
            elif attr in TYPE_CHECK and type(TYPE_CHECK[attr]) != type(value):
                try:
                    value = type(TYPE_CHECK[attr])(value)
                except ValueError:
                    warnings.warn('type check failed for attr {} with value {}'.format(attr, repr(value)))
            elif attr in ['stdev_fragment_size', 'median_fragment_size', 'read_length']:
                try:
                    value = int(value)
                except ValueError:
                    value = float(value)
            elif attr == 'inputs':
                value = value.split(';') if value else []
            section[attr] = value
        config.setdefault(sec, dict()).update(section)

    for lib, section in [(l, config[l]) for l in library_sections]:
        d = dict()
        d.update(defaults)
        d.update(config['qsub'])
        d.update(config['visualization'])
        d.update(config['reference'])
        d.update(section)
        config[lib] = Namespace(**d)
    if len(library_sections) < 1:
        raise UserWarning('configuration file must have 1 or more library sections')
    return [config[l] for l in library_sections]


def main_validate(args):
    """
    - read the evidence
    - assemble contigs from the split reads
    - blat the contigs
    - pair the blatted contigs (where appropriate)
    - TODO: call the breakpoints and summarize the evidence
    """
    FILENAME_PREFIX = re.sub('\.(txt|tsv|tab)$', '', os.path.basename(args.input))
    RAW_EVIDENCE_BAM = os.path.join(args.output, FILENAME_PREFIX + '.raw_evidence.bam')
    CONTIG_BAM = os.path.join(args.output, FILENAME_PREFIX + '.contigs.bam')
    EVIDENCE_BED = os.path.join(args.output, FILENAME_PREFIX + '.evidence.bed')

    PASSED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + PASS_SUFFIX)
    PASSED_BED_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-passed.bed')
    FAILED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-failed.tab')
    CONTIG_BLAT_FA = os.path.join(args.output, FILENAME_PREFIX + '.contigs.fa')
    CONTIG_BLAT_OUTPUT = os.path.join(args.output, FILENAME_PREFIX + '.contigs.blat_out.pslx')
    IGV_BATCH_FILE = os.path.join(args.output, FILENAME_PREFIX + '.igv.batch')
    INPUT_BAM_CACHE = BamCache(args.bam_file, args.stranded_bam)

    evidence_reads = set()  # keep track of collected reads to use for ouput

    split_read_contigs = set()
    chr_to_index = {}
    bpps = read_inputs([args.input])
    evidence_clusters = []
    for bpp in bpps:
        if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            e = GenomeEvidence(
                bpp.break1, bpp.break2,
                INPUT_BAM_CACHE,
                args.reference_genome,
                opposing_strands=bpp.opposing_strands,
                stranded=bpp.stranded,
                untemplated_seq=bpp.untemplated_seq,
                data=bpp.data,
                stdev_fragment_size=args.stdev_fragment_size,
                read_length=args.read_length,
                median_fragment_size=args.median_fragment_size
            )
            evidence_clusters.append(e)
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            e = TranscriptomeEvidence(
                REFERENCE_ANNOTATIONS,
                bpp.break1, bpp.break2,
                INPUT_BAM_CACHE,
                args.reference_genome,
                opposing_strands=bpp.opposing_strands,
                stranded=bpp.stranded,
                untemplated_seq=bpp.untemplated_seq,
                data=bpp.data,
                stdev_fragment_size=args.stdev_fragment_size,
                read_length=args.read_length,
                median_fragment_size=args.median_fragment_size
            )
            evidence_clusters.append(e)
        else:
            raise ValueError('protocol not recognized', bpp.data[COLUMNS.protocol])
    
    for chr, masks in args.masking.items(): # extend masking by read length
        for mask in masks:
            mask.start -= args.read_length
            mask.end += args.read_length
    
    evidence_clusters, filtered_evidence_clusters = filter_on_overlap(evidence_clusters, args.masking)
    
    for i, e in enumerate(evidence_clusters):
        print()
        log(
            '({} of {})'.format(i + 1, len(evidence_clusters)),
            'gathering evidence for:',
            e.data['cluster_id']
        )
        log(e, time_stamp=False)
        log('possible event type(s):', BreakpointPair.classify(e), time_stamp=False)
        log('outer window regions:  {}:{}-{}  {}:{}-{}'.format(
            e.break1.chr, e.outer_window1[0], e.outer_window1[1],
            e.break2.chr, e.outer_window2[0], e.outer_window2[1]), time_stamp=False)
        log('inner window regions:  {}:{}-{}  {}:{}-{}'.format(
            e.break1.chr, e.inner_window1[0], e.inner_window1[1],
            e.break2.chr, e.inner_window2[0], e.inner_window2[1]), time_stamp=False)
        e.load_evidence(log=log)
        log(
            'flanking pairs: {};'.format(len(e.flanking_pairs)),
            'split reads: {}, {};'.format(*[len(a) for a in e.split_reads]),
            'half-mapped reads: {}, {};'.format(*[len(a) for a in e.half_mapped]),
            'spanning-reads: {};'.format(len(e.spanning_reads)),
            'compatible flanking pairs:', len(e.compatible_flanking_pairs),
            time_stamp=False
        )
        e.assemble_contig(log=log)
        log('assembled {} contigs'.format(len(e.contigs)), time_stamp=False)

    log('aligning {} contig sequences'.format(len(blat_sequences)))
    log('will output:', CONTIG_BLAT_FA, CONTIG_BLAT_OUTPUT)
    if len(blat_sequences) > 0:
        blat_contigs(
            evidence,
            INPUT_BAM_CACHE,
            REFERENCE_GENOME=HUMAN_REFERENCE_GENOME,
            blat_2bit_reference=args.blat_2bit_reference,
            blat_fa_input_file=CONTIG_BLAT_FA,
            blat_pslx_output_file=CONTIG_BLAT_OUTPUT,
            clean_files=False
        )
    log('alignment complete')
    event_calls = []
    passes = 0
    
    with open(EVIDENCE_BED, 'w') as fh:
        for index, e in enumerate(evidence):
            fh.write('{}\t{}\t{}\touter-{}\n'.format(
                e.break1.chr, e.outer_window1.start, e.outer_window1.end, e.data[COLUMNS.cluster_id]))
            fh.write('{}\t{}\t{}\touter-{}\n'.format(
                e.break2.chr, e.outer_window2.start, e.outer_window2.end, e.data[COLUMNS.cluster_id]))
            fh.write('{}\t{}\t{}\tinner-{}\n'.format(
                e.break1.chr, e.inner_window1.start, e.inner_window1.end, e.data[COLUMNS.cluster_id]))
            fh.write('{}\t{}\t{}\tinner-{}\n'.format(
                e.break2.chr, e.inner_window2.start, e.inner_window2.end, e.data[COLUMNS.cluster_id]))
            print()
            log('({} of {}) calling events for:'.format
                (index + 1, len(evidence)), e.data[COLUMNS.cluster_id], e.putative_event_types())
            log('source:', e, time_stamp=False)
            calls = []
            failure_comment = None
            try:
                calls.extend(call_events(e))
                event_calls.extend(calls)
            except UserWarning as err:
                log('warning: error in calling events', repr(err), time_stamp=False)
                failure_comment = str(err)
            if len(calls) == 0:
                failure_comment = ['zero events were called'] if failure_comment is None else failure_comment
                e.data[COLUMNS.filter_comment] = failure_comment
                filtered_evidence_clusters.append(e)
            else:
                passes += 1

            log('called {} event(s)'.format(len(calls)))
            for ev in calls:
                log(ev, time_stamp=False)
                log(ev.event_type, ev.call_method, time_stamp=False)
                log('remapped reads: {}; split reads: [{}, {}], flanking pairs: {}'.format(
                    0 if not ev.contig else len(ev.contig.input_reads),
                    len(ev.break1_split_reads), len(ev.break2_split_reads),
                    len(ev.flanking_pairs)), time_stamp=False)

    if len(failed_cluster_rows) + passes != len(evidence):
        raise AssertionError(
            'totals do not match pass + fails == total', passes, len(failed_cluster_rows), len(evidence))
    # write the output validated clusters (split by type and contig)
    validation_batch_id = build_batch_id(prefix='validation-')
    for i, ec in enumerate(event_calls):
        b1_homseq = None
        b2_homseq = None
        try:
            b1_homseq, b2_homseq = ec.breakpoint_sequence_homology(HUMAN_REFERENCE_GENOME)
        except AttributeError:
            pass
        ec.data.update({
            COLUMNS.validation_id: '{}-{}'.format(validation_batch_id, i + 1),
            COLUMNS.break1_homologous_seq: b1_homseq,
            COLUMNS.break2_homologous_seq: b2_homseq,
        })

    output_tabbed_file(event_calls, PASSED_OUTPUT_FILE)
    output_tabbed_file(filtered_evidence_clusters, FAILED_OUTPUT_FILE)

    with pysam.AlignmentFile(CONTIG_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        log('writing:', CONTIG_BAM)
        for ev in evidence:
            for c in ev.contigs:
                for read1, read2 in c.alignments:
                    read1.cigar = cigar_tools.convert_for_igv(read1.cigar)
                    fh.write(read1)
                    if read2:
                        read2.cigar = cigar_tools.convert_for_igv(read2.cigar)
                        fh.write(read2)

    # write the evidence
    with pysam.AlignmentFile(EVIDENCE_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        log('writing:', EVIDENCE_BAM)
        reads = set()
        for ev in evidence:
            temp = ev.supporting_reads()
            reads.update(temp)
        for read in reads:
            read.cigar = cigar_tools.convert_for_igv(read.cigar)
            fh.write(read)
    # now sort the contig bam
    sort = re.sub('.bam$', '.sorted.bam', CONTIG_BAM)
    log('sorting the bam file:', CONTIG_BAM)
    if SAMTOOLS_VERSION[0] < 1:
        subprocess.call(samtools_v0_sort(CONTIG_BAM, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(CONTIG_BAM, sort), shell=True)
    CONTIG_BAM = sort
    log('indexing the sorted bam:', CONTIG_BAM)
    subprocess.call(['samtools', 'index', CONTIG_BAM])

    # then sort the evidence bam file
    sort = re.sub('.bam$', '.sorted.bam', EVIDENCE_BAM)
    log('sorting the bam file:', EVIDENCE_BAM)
    if SAMTOOLS_VERSION[0] < 1:
        subprocess.call(samtools_v0_sort(EVIDENCE_BAM, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(EVIDENCE_BAM, sort), shell=True)
    EVIDENCE_BAM = sort
    log('indexing the sorted bam:', EVIDENCE_BAM)
    subprocess.call(['samtools', 'index', EVIDENCE_BAM])

    # write the igv batch file
    with open(IGV_BATCH_FILE, 'w') as fh:
        log('writing:', IGV_BATCH_FILE)
        fh.write('new\ngenome {}\n'.format(args.reference_genome))

        fh.write('load {} name="{}"\n'.format(PASSED_BED_FILE, 'passed events'))
        fh.write('load {} name="{}"\n'.format(CONTIG_BAM, 'aligned contigs'))
        fh.write('load {} name="{}"\n'.format(EVIDENCE_BED, 'evidence windows'))
        fh.write('load {} name="{}"\n'.format(EVIDENCE_BAM, 'raw evidence'))
        fh.write('load {} name="{} {} input"\n'.format(args.bam_file, args.library, args.protocol))


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
    parser.add_argument('output', help='path to the output directory')
    parser.add_argument(
        '-f', '--force_overwrite', default=False, type=TSV.tsv_boolean,
        help='set flag to overwrite existing reviewed files'
    )
    
    if pstep == PIPELINE_STEP.ANNOTATE:
        parser.add_argument(
            '--output_svgs', default=True, type=TSV.tsv_boolean,
            help='set flag to suppress svg drawings of putative annotations')
        parser.add_argument(
            '--min_orf_size', default=MIN_ORF_SIZE, type=int, help='minimum size for putative ORFs'
        )
        parser.add_argument(
            '--max_orf_cap', default=MAX_ORF_CAP, type=int, help='keep the n longest orfs'
        )
        parser.add_argument(
            '--min_domain_mapping_match', default=0.8, type=float,
            help='minimum percent match for the domain to be considered aligned'
        )
        g = parser.add_argument_group('visualization options')
        g.add_argument(
            '--domain_regex_filter', default='^PF\d+$',
            help='only show domains which names (external identifiers) match the given pattern'
        )

    if pstep == PIPELINE_STEP.CLUSTER or pstep == PIPELINE_STEP.ANNOTATE:
        parser.add_argument(
            '--max_proximity', default=MAX_PROXIMITY, type=int,
            help='The maximum distance away from breakpoints to look for proximal genes'
        )
        parser.add_argument('-n', '--inputs', required=True, nargs='+', help='1 or more input files')
    elif pstep == PIPELINE_STEP.VALIDATE or pstep == PIPELINE_STEP.PAIR:
        parser.add_argument('-n', '--input', help='path to the input file', required=True)
    
    if pstep == PIPELINE_STEP.CLUSTER or pstep == PIPELINE_STEP.VALIDATE:
        parser.add_argument('library', help='library name')
        parser.add_argument('protocol', help='the library protocol: genome or transcriptome', choices=PROTOCOL.values())
    
    if pstep == PIPELINE_STEP.CLUSTER:
        g = parser.add_argument_group('output arguments')
        g.add_argument(
            '--max_files', default=MAX_FILES, type=int, dest='max_files',
            help='defines the maximum number of files that can be created')
        g.add_argument(
            '--min_clusters_per_file', default=MIN_CLUSTERS_PER_FILE, type=int,
            help='defines the minimum number of clusters per file')
        parser.add_argument(
            '-r', '--cluster_radius', help='radius to use in clustering', default=CLUSTER_RADIUS, type=int)
        parser.add_argument(
            '-k', '--cluster_clique_size',
            help='parameter used for computing cliques, smaller is faster, above 20 will be slow',
            default=CLUSTER_CLIQUE_SIZE, type=int
        )
        parser.add_argument(
            '--uninformative_filter', default=True, help='If flag is False then the clusters will not be filtered '
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

        parser.add_argument('-p', '--protocol', required=True, choices=PROTOCOL.values())

    args = parser.parse_args()
    if pstep == PIPELINE_STEP.VALIDATE:
        args.samtools_version = get_samtools_version()
        args.blat_version = get_blat_version()
    args.output = os.path.abspath(args.output)
    if os.path.exists(args.output) and not args.force_overwrite:
        parser.print_help()
        print(
            '\nerror: output directory {} exists, --force_overwrite must be specified or the directory removed'.format(
                repr(args.output)))
        exit(1)
    return parser, args


def filter_on_overlap(bpps, regions_by_reference_name):
    log('filtering on overlaps with regions')
    failed = []
    passed = []
    for bpp in bpps:
        overlaps = False
        for r in regions_by_reference_name.get(bpp.break1.chr, []):
            if Interval.overlaps(r, bpp.break1):
                overlaps = True
                bpp.data['failure_comment'] = 'overlapped masked region: ' + str(r)
                break
        for r in regions_by_reference_name.get(bpp.break2.chr, []):
            if overlaps:
                break
            if Interval.overlaps(r, bpp.break2):
                overlaps = True
                bpp.data[COLUMNS.filter_comment] = 'overlapped masked region: ' + str(r)
        if overlaps:
            failed.append(bpp)
        else:
            passed.append(bpp)
    log('filtered', len(bpps), 'to', len(passed))
    return passed, failed


def main_pipeline(parser, args):

    READ_FILES = {}

    if os.path.exists(args.output) and not args.force_overwrite:
        print('error: must specify the overwrite option or delete the existing output directory')
        parser.print_help()
        exit(1)
    args.output = os.path.abspath(args.output)
    configs = read_config(args.config)
    # read the config
    # set up the directory structure and run svmerge
    annotation_jobs = []
    for sec in configs:
        base = os.path.join(args.output, '{}_{}'.format(sec.library, sec.protocol))
        log('setting up the directory structure for', sec.library, 'as', base)
        base = os.path.join(args.output, '{}_{}'.format(sec.library, sec.protocol))
        cluster_output = mkdirp(os.path.join(base, 'clustering'))
        validation_output = mkdirp(os.path.join(base, 'validation'))
        annotation_output = mkdirp(os.path.join(base, 'annotation'))

        # run the merge
        log('clustering')
        setattr(sec, 'output', os.path.join(base, 'clustering'))
        merge_args = {
            'output': cluster_output,
            'max_proximity': sec.max_proximity,
            'cluster_radius': CLUSTER_RADIUS,
            'cluster_clique_size': CLUSTER_CLIQUE_SIZE,
            'max_files': MAX_FILES,
            'min_clusters_per_file': MIN_CLUSTERS_PER_FILE,
            'uninformative_filter': True
        }
        merge_args.update(sec.__dict__)
        ann = merge_args['annotations']
        merge_args['annotations'] = READ_FILES.get(ann, ann)
        merge_args = Namespace(**merge_args)
        output_files = cluster_main(merge_args)
        READ_FILES[ann] = getattr(merge_args, 'annotations')
        merge_file_prefix = None
        for f in output_files:
            m = re.match('^(?P<prefix>.*\D)\d+.tab$', f)
            if not m:
                raise UserWarning('cluster file did not match expected format', f)
            if merge_file_prefix is None:
                merge_file_prefix = m.group('prefix')
            elif merge_file_prefix != m.group('prefix'):
                raise UserWarning('merge file prefixes are not consistent', output_files)

        # now set up the qsub script for the validation and the held job for the annotation
        validation_args = {
            'output': validation_output,
            'masking': sec.masking,
            'reference_genome': sec.reference_genome,
            'blat_2bit_reference': sec.blat_2bit_reference,
            'annotations': sec.annotations,
            'library': sec.library,
            'bam_file': sec.bam_file,
            'protocol': sec.protocol,
            'read_length': sec.read_length,
            'stdev_fragment_size': sec.stdev_fragment_size,
            'median_fragment_size': sec.median_fragment_size,
            'force_overwrite': args.force_overwrite,
            'stranded_bam': sec.stranded_bam
        }
        for attr in sorted(VALIDATION_DEFAULTS.__dict__.keys()):
            validation_args[attr] = getattr(sec, attr)

        qsub = os.path.join(validation_output, 'qsub.sh')
        validation_jobname = 'validation_{}_{}'.format(sec.library, sec.protocol)
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=sec.queue, memory=sec.memory, name=validation_jobname, output=validation_output
                ) + '\n')
            fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            temp = ['--{} {}'.format(k, v) for k, v in validation_args.items() if not isinstance(v, str) and v is not None]
            temp.extend(['--{} "{}"'.format(k, v) for k, v in validation_args.items() if isinstance(v, str) and v is not None])
            validation_args = temp
            validation_args.append('-n {}$SGE_TASK_ID.tab'.format(merge_file_prefix))
            fh.write('python {}/mavis_validate.py {}\n'.format(basedir, ' \\\n\t'.join(validation_args)))

        # set up the annotations job
        # for all files with the right suffix
        annotation_args = {
            'output': annotation_output,
            'reference_genome': sec.reference_genome,
            'annotations': sec.annotations,
            'template_metadata': sec.template_metadata,
            'min_orf_size': sec.min_orf_size,
            'max_orf_cap': sec.max_orf_cap,
            'min_domain_mapping_match': sec.min_domain_mapping_match,
            'domain_regex_filter': sec.domain_regex_filter,
            'max_proximity': sec.max_proximity
        }
        temp = ['--{} {}'.format(k, v) for k, v in annotation_args.items() if not isinstance(v, str) and v is not None]
        temp.extend(['--{} "{}"'.format(k, v) for k, v in annotation_args.items() if isinstance(v, str) and v is not None])
        annotation_args = temp
        annotation_args.append('--input {}/*{}'.format(validation_output, PASS_SUFFIX))
        qsub = os.path.join(annotation_output, 'qsub.sh')
        annotation_jobname = 'annotation_{}_{}'.format(sec.library, sec.protocol)
        annotation_jobs.append(annotation_jobname)
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=sec.queue, memory=sec.memory, name=annotation_jobname, output=annotation_output
                ) + '\n')
            fh.write('#$ -hold_jid {}\n'.format(validation_jobname))
            fh.write('python {}/mavis_annotate.py {}\n'.format(basedir, ' \\\n\t'.join(annotation_args)))

    # set up scripts for the pairing held on all of the annotation jobs
    pairing_output = mkdirp(os.path.join(base, 'pairing'))
    qsub = os.path.join(pairing_output, 'qsub.sh')
    with open(qsub, 'w') as fh:
        log('writing:', qsub)


def read_inputs(inputs, force_stranded=False, **kwargs):
    bpps = []
    kwargs.setdefault('require', [])
    kwargs['require'] = list(set(kwargs['require'] + [COLUMNS.library, COLUMNS.protocol]))
    kwargs.setdefault('validate', {})
    kwargs['validate'][COLUMNS.library] = '^[\w-]+$'
    kwargs.setdefault('_in', {})
    kwargs['_in'][COLUMNS.protocol] = PROTOCOL
    for finput in inputs:
        log('loading:', finput)
        bpps.extend(read_bpp_from_input_file(
            finput, force_stranded=force_stranded,
            **kwargs
        ))
    return bpps


def output_tabbed_file(bpps, filename):
    header = set()
    rows = []
    for row in bpps:
        try:
            row = row.flatten()
        except AttributeError:
            pass
        rows.append(row)
        header.update(row)

    header = sort_columns(header)
    
    with open(filename, 'w') as fh:
        fh.write('#' + '\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row[c]) for c in header]) + '\n')


def main_cluster(args):
    # output files
    cluster_batch_id = build_batch_id(prefix='cluster-')
    UNINFORM_OUTPUT = os.path.join(args.output, 'uninformative_clusters.txt')
    CLUSTER_ASSIGN_OUTPUT = os.path.join(args.output, 'cluster_assignment.tab')
    CLUSTER_ASSIGN_OUTPUT = os.path.join(args.output, 'cluster_assignment.tab')
    CLUSTER_BED_OUTPUT = os.path.join(args.output, 'clusters.bed')
    split_file_name_func = lambda x: os.path.join(args.output, '{}-{}.tab'.format(cluster_batch_id, x))
    # load the input files
    breakpoint_pairs = read_inputs(
        args.inputs, args.stranded_bam,
        cast={COLUMNS.tools: lambda x: set(';'.split(x)) if x else set()},
        add={COLUMNS.library: args.library, COLUMNS.protocol: args.protocol}
    )
    # filter by masking file
    breakpoint_pairs, filtered_bpp = filter_on_overlap(breakpoint_pairs, args.masking)

    log('computing clusters')
    clusters = cluster_breakpoint_pairs(breakpoint_pairs, r=args.cluster_radius, k=args.cluster_clique_size)

    hist = {}
    length_hist = {}
    for cluster, input_pairs in clusters.items():
        hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
        c1 = round(len(cluster[0]), -2)
        c2 = round(len(cluster[1]), -2)
        length_hist[c1] = length_hist.get(c1, 0) + 1
        length_hist[c2] = length_hist.get(c2, 0) + 1
        cluster.data[COLUMNS.cluster_id] = '{}-{}'.format(cluster_id_prefix, cluster_id)
        cluster.data[COLUMNS.cluster_size] = len(input_pairs)
        temp = set()
        for p in input_pairs:
            temp.update(p.data[COLUMNS.tools])
        cluster.data[COLUMNS.tools] = ';'.join(sorted(list(temp)))
        cluster_id += 1
    log('computed', len(clusters), 'clusters', time_stamp=False)
    log('cluster input pairs distribution', sorted(hist.items()), time_stamp=False)
    log('cluster intervals lengths', sorted(length_hist.items()), time_stamp=False)
    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    mkdirp(args.output)
    
    with open(CLUSTER_ASSIGN_OUTPUT, 'w') as fh:
        header = set()
        log('writing:', CLUSTER_ASSIGN_OUTPUT)
        rows = {}

        for cluster, input_pairs in clusters.items():
            for p in input_pairs:
                if p in rows:
                    rows[p][COLUMNS.tools].update(p.data[COLUMNS.tools])
                else:
                    rows[p] = BreakpointPair.flatten(p)
                rows[p].setdefault('clusters', set()).add(cluster.data[COLUMNS.cluster_id])
        for row in rows.values():
            row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
            row[COLUMNS.tools] = ';'.join(sorted(list(row[COLUMNS.tools])))
            row[COLUMNS.library] = args.library
            row[COLUMNS.protocol] = args.protocol
        output_tabbed_file(rows, CLUSTER_ASSIGN_OUTPUT)

    output_files = []
    # filter clusters based on annotations
    # decide on the number of clusters to validate per job
    pass_clusters = clusters
    fail_clusters = []

    if args.uninformative_filter:
        for cluster in clusters:
            # loop over the annotations
            overlaps_gene = False
            w1 = Interval(cluster.break1.start - args.max_proximity, cluster.break1.end + args.max_proximity)
            w2 = Interval(cluster.break2.start - args.max_proximity, cluster.break2.end + args.max_proximity)
            for gene in args.annotations.get(cluster.break1.chr, []):
                if Interval.overlaps(gene, w1):
                    overlaps_gene = True
                    break
            for gene in args.annotations.get(cluster.break2.chr, []):
                if Interval.overlaps(gene, w2):
                    overlaps_gene = True
                    break
            if overlaps_gene:
                pass_clusters.append(cluster)
            else:
                fail_clusters.append(cluster)
    assert(len(fail_clusters) + len(pass_clusters) == len(clusters))
    log('filtered', len(fail_clusters), 'clusters as not informative')
    output_tabbed_file(fail_clusters, UNINFORM_OUTPUT)

    JOB_SIZE = args.min_clusters_per_file
    if len(pass_clusters) // args.min_clusters_per_file > args.max_files - 1:
        JOB_SIZE = len(pass_clusters) // args.max_files
        assert(len(pass_clusters) // JOB_SIZE == args.max_files)
    
    bedfile = os.path.join(args.output, 'clusters.bed')
    log('writing:', bedfile)
    write_bed_file(bedfile, clusters)

    job_ranges = list(range(0, len(pass_clusters), JOB_SIZE))
    if job_ranges[-1] != len(pass_clusters):
        job_ranges.append(len(pass_clusters))
    job_ranges = zip(job_ranges, job_ranges[1::])

    for i, jrange in enumerate(job_ranges):
        # generate an output file
        filename = split_file_name_func(i + 1)
        output_files.append(filename)
        output_tabbed_file(pass_clusters[jrange[0]:jrange[1]], filename)

    return output_files

def main_annotate(args):
    DRAWINGS_DIRECTORY = os.path.join(args.output, 'drawings')
    TABBED_OUTPUT_FILE = os.path.join(args.output, 'annotations.tab')
    FA_OUTPUT_FILE = os.path.join(args.output, 'annotations.fusion-cdna.fa')
    
    mkdirp(DRAWINGS_DIRECTORY)
    # test that the sequence makes sense for a random transcript
    bpps = read_inputs(args.inputs, _in={COLUMNS.protocol: PROTOCOL})
    log('read {} breakpoint pairs'.format(len(bpps)))

    annotations = annotate_events(
        bpps, 
        reference_genome=args.reference_genome, annotations=args.annotations, 
        min_orf_size=args.min_orf_size, min_domain_mapping_match=args.min_domain_mapping_match,
        max_orf_cap=args.max_orf_cap
    ) 

    id_prefix = build_batch_id(prefix='annotation', suffix='-')
    rows = []  # hold the row information for the final tsv file
    fa_sequences = {}
    # now try generating the svg
    DS = DiagramSettings()
    DS.DOMAIN_NAME_REGEX_FILTER = args.domain_regex_filter

    for i, ann in enumerate(annotations):
        ann.data[COLUMNS.annotation_id] = id_prefix + str(i + 1)
        row = ann.flatten()
        row[COLUMNS.break1_strand] = ann.transcript1.get_strand()
        row[COLUMNS.break2_strand] = ann.transcript2.get_strand()
        row[COLUMNS.fusion_sequence_fasta_file] = FA_OUTPUT_FILE

        log('current annotation', annotation_id, ann.transcript1, ann.transcript2, ann.event_type)

        # try building the fusion product
        ann_rows = []
        # add fusion information to the current row
        for t in ft.transcripts:
            fusion_fa_id = '{}_{}'.format(annotation_id, t.splicing_pattern.splice_type)
            if fusion_fa_id in fa_sequences:
                raise AssertionError('should not be duplicate fa sequence ids', fusion_fa_id)
            fa_sequences[fusion_fa_id] = ft.get_cdna_seq(t.splicing_pattern)

            # duplicate the row for each translation
            for tl in t.translations:
                nrow = dict()
                nrow.update(row)
                nrow[COLUMNS.fusion_splicing_pattern] = tl.transcript.splicing_pattern.splice_type
                nrow[COLUMNS.fusion_cdna_coding_start] = tl.start
                nrow[COLUMNS.fusion_cdna_coding_end] = tl.end
                nrow[COLUMNS.fusion_sequence_fasta_id] = fusion_fa_id

                domains = []
                for dom in tl.domains:
                    m, t = dom.score_region_mapping()
                    temp = {
                        "name": dom.name,
                        "sequences": dom.get_seqs(),
                        "regions": [
                            {"start": dr.start, "end": dr.end} for dr in sorted(dom.regions, key=lambda x: x.start)
                        ],
                        "mapping_quality": round(m * 100 / t, 0),
                        "matches": m
                    }
                    domains.append(temp)
                nrow[COLUMNS.fusion_mapped_domains] = json.dumps(domains)
                ann_rows.append(nrow)

        drawing = None
        retry_count = 0
        while drawing is None:  # continue if drawing error and increase width
            try:
                canvas, legend = draw_sv_summary_diagram(
                    DS, ann, REFERENCE_GENOME=REFERENCE_GENOME, templates=TEMPLATES)

                gene_aliases1 = 'NA'
                gene_aliases2 = 'NA'
                try:
                    if len(ann.transcript1.gene.aliases) > 0:
                        gene_aliases1 = '-'.join(ann.transcript1.gene.aliases)
                    if ann.transcript1.is_best_transcript:
                        gene_aliases1 = 'b-' + gene_aliases1
                except AttributeError:
                    pass
                try:
                    if len(ann.transcript2.gene.aliases) > 0:
                        gene_aliases2 = '-'.join(ann.transcript2.gene.aliases)
                    if ann.transcript2.is_best_transcript:
                        gene_aliases2 = 'b-' + gene_aliases2
                except AttributeError:
                    pass
                try:
                    if determine_prime(ann.transcript1, ann.break1) == PRIME.THREE:
                        gene_aliases1, gene_aliases2 = gene_aliases2, gene_aliases1
                except NotSpecifiedError:
                    pass

                name = 'chr{}_chr{}.{}_{}.{}'.format(
                    ann.break1.chr, ann.break2.chr, gene_aliases1, gene_aliases2, ann.data[COLUMNS.annotation_id])

                drawing = os.path.join(DRAWINGS_DIRECTORY, name + '.svg')
                l = os.path.join(DRAWINGS_DIRECTORY, name + '.legend.json')
                for r in ann_rows + [row]:
                    r[COLUMNS.annotation_figure] = drawing
                    r[COLUMNS.annotation_figure_legend] = l
                log('generating svg:', drawing, time_stamp=False)
                canvas.saveas(drawing)

                log('generating legend:', l, time_stamp=False)
                with open(l, 'w') as fh:
                    json.dump(legend, fh)
                break
            except DrawingFitError as err:
                log('extending width:', DS.WIDTH, DS.WIDTH + 500, time_stamp=False)
                DS.WIDTH += 500
                retry_count += 1
                if retry_count > 10:
                    raise err
        if len(ann_rows) == 0:
            rows.append(row)
        else:
            rows.extend(ann_rows)
    
    output_tabbed_file(rows, TABBED_OUTPUT_FILE)

    with open(FA_OUTPUT_FILE, 'w') as fh:
        log('writing:', FA_OUTPUT_FILE)
        for name, seq in sorted(fa_sequences.items()):
            fh.write('> {}\n'.format(name))
            fh.write('{}\n'.format(seq))

def main_pairing(args):
    pass

def main_summary(args):
    pass

def main():
    def usage(err):
        name = os.path.basename(__file__)
        print('usage: {} {{cluster,validate,annotate,pairing,summary,pipeline}}'.format(name))
        print('{}: error:'.format(name), err)
        exit(1)
    
    if len(sys.argv) < 2:
        usage('the <pipeline step> argument is required')
    
    pstep = sys.argv.pop(1)
    sys.argv[0] = '{} {}'.format(sys.argv[0], pstep)
    
    if pstep == PIPELINE_STEP.SUMMARY:
        raise NotImplementedError('summary script has not been written')
    try:
        args = parse_arguments(pstep)
        config = None
        log('input arguments')
        for arg, val in sorted(args.__dict__.items()):
            log(arg, '=', repr(val), time_stamp=False)
        if pstep == PIPELINE_STEP.PIPELINE:
            config = read_config(args.config)
            args = config[0]
        
        # load the reference files if they have been given?
        if any([
                pstep not in [PIPELINE_STEP.PIPELINE, PIPELINE.VALIDATE],
                args.uninformative_filter,
                pstep == PIPELINE_STEP.VALIDATE and args.protocol == PROTOCOL.TRANS]):
            log('loading:', args.annotations)
            args.__dict__['annotations'] = load_reference_genes(args.annotations)
        try:
            log('loading:' if not args.low_memory else 'indexing:', args.reference_genome)
            args.__dict__['reference_genome'] = load_reference_genome(args.reference_genome, args.low_memory)
        except AttributeError:
            pass
        try:
            log('loading:', args.masking)
            args.__dict__['masking'] = load_masking_regions(args.masking)
        except AttributeError:
            pass
        try:
            log('loading:', args.template_metadata)
            args.__dict__['template_metadata'] = load_templates(args.template_metadata)
        except AttributeError:
            pass
        
        if pstep == PIPELINE_STEP.CLUSTER:
            main_cluster(args)
        elif pstep == PIPELINE_STEP.VALIDATE:
            main_validate(args)
        elif pstep == PIPELINE_STEP.ANNOTATE:
            main_annotate(args)
        elif pstep == PIPELINE_STEP.PAIR:
            main_pairing(args)
        elif pstep == PIPELINE_STEP.SUMMARY:
            main_summary(args)
        else:  # PIPELINE
            main_pipeline(config)


    except KeyError:
        usage('unrecognized value {} for <pipeline step>'.format(repr(pstep)))

if __name__ == '__main__':
    main()
