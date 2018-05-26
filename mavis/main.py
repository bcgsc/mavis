#!python
import argparse
import platform
import os
import time

import tab

from . import __version__
from .align import get_aligner_version
from . import annotate as _annotate
from .annotate import main as annotate_main
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .cluster import main as cluster_main
from .config import augment_parser, MavisConfig, generate_config, get_metavar, CustomHelpFormatter, RangeAppendAction
from .constants import SUBCOMMAND, PROTOCOL, float_fraction
from .error import DrawingFitError
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS, DiagramSettings
from .illustrate.diagram import draw_multi_transcript_overlay
from .illustrate.scatter import bam_to_scatter
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .pairing import main as pairing_main
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .summary import main as summary_main
from .tools import convert_tool_output, SUPPORTED_TOOL
from .util import bash_expands, get_env_variable, LOG, log_arguments, MavisNamespace, mkdirp, output_tabbed_file
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from .validate import main as validate_main
from . import schedule as _schedule


PROGNAME = 'mavis'
EXIT_OK = 0
EXIT_ERROR = 1


def check_overlay_args(args, parser):
    """
    parse the overlay options and check the formatting
    """
    # check complex options
    for marker in args.markers:
        if len(marker) < 3:
            marker.append(marker[-1])
        try:
            marker[1] = int(marker[1])
            marker[2] = int(marker[2])
        except ValueError:
            parser.error('argument --marker: start and end must be integers: {}'.format(marker))

    defaults = [None, None, 0.5, None, True]
    bam_file, density, ymax, stranded = range(1, 5)

    for plot in args.read_depth_plots:
        for i, d in enumerate(defaults):
            if i >= len(plot):
                plot.append(d)
        if not os.path.exists(plot[bam_file]):
            parser.error('argument --read_depth_plots: the bam file given does not exist: {}'.format(plot[bam_file]))
        try:
            plot[density] = float(plot[density])
            if plot[density] < 0 or plot[density] > 1:
                raise ValueError()
        except ValueError:
            parser.error('argument --read_depth_plots: density must be an float between 0 and 1: {}'.format(plot[density]))
        try:
            if str(plot[ymax]).lower() in ['null', 'none']:
                plot[ymax] = None
            else:
                plot[ymax] = int(plot[ymax])
        except ValueError:
            parser.error('argument --read_depth_plots: ymax must be an integer: {}'.format(plot[ymax]))
        try:
            plot[stranded] = tab.cast_boolean(plot[stranded])
        except TypeError:
            parser.error('argument --read_depth_plots: stranded must be an boolean: {}'.format(plot[stranded]))
    return args


def overlay_main(
    gene_name, output, buffer_length, read_depth_plots, markers,
    annotations,
    drawing_width_iter_increase, max_drawing_retries, min_mapping_quality,
    ymax_color='#FF0000',
    **kwargs
):
    """
    generates an overlay diagram
    """
    annotations.load()
    # check options formatting
    gene_to_draw = None

    for chrom in annotations.content:
        for gene in annotations.content[chrom]:
            if gene_name in gene.aliases or gene_name == gene.name:
                gene_to_draw = gene
                LOG('Found target gene: {}(aka. {}) {}:{}-{}'.format(gene.name, gene.aliases, gene.chr, gene.start, gene.end))
                break
    if gene_to_draw is None:
        raise KeyError('Could not find gene alias or id in annotations file', gene_name)

    settings = DiagramSettings(**kwargs)

    genomic_min = max(gene_to_draw.start - buffer_length, 1)
    genomic_max = gene_to_draw.end + buffer_length

    plots = []
    for axis_name, bam_file, density, ymax, stranded in read_depth_plots:
        # one plot per bam
        plots.append(bam_to_scatter(
            bam_file, gene_to_draw.chr, genomic_min, genomic_max,
            strand=gene_to_draw.get_strand() if stranded else None,
            ymax=ymax,
            density=density,
            axis_name=axis_name,
            min_mapping_quality=min_mapping_quality,
            ymax_color=ymax_color
        ))

    for i, (marker_name, marker_start, marker_end) in enumerate(markers):
        markers[i] = _annotate.base.BioInterval(gene_to_draw.chr, marker_start, marker_end, name=marker_name)

    canvas = None
    attempts = 1
    while True:
        try:
            canvas = draw_multi_transcript_overlay(
                settings, gene_to_draw,
                vmarkers=markers,
                plots=plots,
                window_buffer=buffer_length,
                log=LOG
            )
            break
        except DrawingFitError as err:
            if attempts > max_drawing_retries:
                raise err
            LOG('Drawing fit: extending window', drawing_width_iter_increase)
            settings.width += drawing_width_iter_increase
            attempts += 1

    svg_output_file = os.path.join(output, '{}_{}_overlay.svg'.format(gene_to_draw.name, gene_name))
    LOG('writing:', svg_output_file)

    canvas.saveas(svg_output_file)


def convert_main(inputs, outputfile, file_type, strand_specific=False, assume_no_untemplated=True):
    bpp_results = convert_tool_output(inputs, file_type, strand_specific, LOG, True, assume_no_untemplated=assume_no_untemplated)
    if os.path.dirname(outputfile):
        mkdirp(os.path.dirname(outputfile))
    output_tabbed_file(bpp_results, outputfile)



def main():
    """
    sets up the parser and checks the validity of command line args
    loads reference files and redirects into subcommand main functions
    """
    start_time = int(time.time())

    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    augment_parser(['version'], parser)
    subp = parser.add_subparsers(dest='command', help='specifies which step/stage in the pipeline or which subprogram to use')
    subp.required = True
    required = {}  # hold required argument group by subparser command name
    optional = {}  # hold optional argument group by subparser command name
    for command in SUBCOMMAND.values():
        subparser = subp.add_parser(command, formatter_class=CustomHelpFormatter, add_help=False)
        required[command] = subparser.add_argument_group('required arguments')
        optional[command] = subparser.add_argument_group('optional arguments')
        augment_parser(['help', 'version'], optional[command])

    # config arguments
    required[SUBCOMMAND.CONFIG].add_argument('-w', '--write', help='path to the new configuration file', required=True, metavar='FILEPATH')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--library',
        metavar='<name> {genome,transcriptome} {diseased,normal} [strand_specific] [/path/to/bam/file]',
        action=RangeAppendAction, help='configuration for libraries to be analyzed by mavis', nmin=3, nmax=5
    )
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--input', help='path to an input file or filter for mavis followed by the library names it '
        'should be used for', nmin=2, action=RangeAppendAction, metavar='FILEPATH <name> [<name> ...]'
    )
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--assign', help='library name followed by path(s) to input file(s) or filter names. This represents the list'
        ' of inputs that should be used for the library', action=RangeAppendAction, nmin=2,
        metavar='<name> FILEPATH [FILEPATH ...]')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--genome_bins', default=get_env_variable('genome_bins', 100), type=int, metavar=get_metavar(int),
        help='number of bins/samples to use in calculating the fragment size stats for genomes')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--transcriptome_bins', default=get_env_variable('transcriptome_bins', 500), type=int, metavar=get_metavar(int),
        help='number of genes to use in calculating the fragment size stats for genomes')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--distribution_fraction', default=get_env_variable('distribution_fraction', 0.97), type=float_fraction, metavar=get_metavar(float),
        help='the proportion of the distribution of calculated fragment sizes to use in determining the stdev')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--convert', nmin=3,
        metavar='<alias> FILEPATH [FILEPATH ...] {{{}}} [stranded]'.format(','.join(SUPPORTED_TOOL.values())),
        help='input file conversion for internally supported tools', action=RangeAppendAction)
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--external_conversion', metavar=('<alias>', '<"command">'), nargs=2, default=[],
        help='alias for use in inputs and full command (quoted)', action='append')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--add_defaults', default=False, action='store_true', help='write current defaults for all non-specified options to the config output')
    augment_parser(['annotations'], optional[SUBCOMMAND.CONFIG])
    # add the optional annotations file (only need this is auto generating bam stats for the transcriptome)
    augment_parser(['skip_stage'], optional[SUBCOMMAND.CONFIG])

    # convert
    required[SUBCOMMAND.CONVERT].add_argument(
        '--file_type', choices=sorted(SUPPORTED_TOOL.values()),
        required=True, help='Indicates the input file type to be parsed')
    augment_parser(['strand_specific', 'assume_no_untemplated'], optional[SUBCOMMAND.CONVERT])
    required[SUBCOMMAND.CONVERT].add_argument('--outputfile', '-o', required=True, help='path to the outputfile', metavar='FILEPATH')

    for command in set(SUBCOMMAND.values()) - {SUBCOMMAND.CONFIG, SUBCOMMAND.CONVERT}:
        required[command].add_argument('-o', '--output', help='path to the output directory', required=True)

    # pipeline
    augment_parser(['config'], required[SUBCOMMAND.PIPELINE])
    optional[SUBCOMMAND.PIPELINE].add_argument(
        '--skip_stage', choices=[SUBCOMMAND.CLUSTER, SUBCOMMAND.VALIDATE], action='append', default=[],
        help='Use flag once per stage to skip. Can skip clustering or validation or both')
    optional[SUBCOMMAND.SCHEDULE].add_argument('--submit', action='store_true', default=False, help='submit jobs to the the scheduler specified')

    # add the inputs argument
    for command in [SUBCOMMAND.CLUSTER, SUBCOMMAND.ANNOTATE, SUBCOMMAND.VALIDATE, SUBCOMMAND.PAIR, SUBCOMMAND.SUMMARY, SUBCOMMAND.CONVERT]:
        required[command].add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')

    # cluster
    augment_parser(
        ['library', 'protocol', 'strand_specific', 'disease_status'],
        required[SUBCOMMAND.CLUSTER])
    augment_parser(list(CLUSTER_DEFAULTS.keys()) + ['masking', 'annotations'], optional[SUBCOMMAND.CLUSTER])

    # validate
    augment_parser(
        [
            'library', 'protocol', 'bam_file', 'read_length', 'stdev_fragment_size',
            'median_fragment_size', 'strand_specific', 'reference_genome', 'aligner_reference'
        ],
        required[SUBCOMMAND.VALIDATE]
    )
    augment_parser(VALIDATION_DEFAULTS.keys(), optional[SUBCOMMAND.VALIDATE])
    augment_parser(['masking', 'annotations'], optional[SUBCOMMAND.VALIDATE])

    # annotate
    augment_parser(
        ['library', 'protocol', 'annotations', 'reference_genome'],
        required[SUBCOMMAND.ANNOTATE]
    )
    augment_parser(['max_proximity', 'masking', 'template_metadata'], optional[SUBCOMMAND.ANNOTATE])
    augment_parser(list(_annotate.constants.DEFAULTS.keys()) + list(ILLUSTRATION_DEFAULTS.keys()), optional[SUBCOMMAND.ANNOTATE])

    # pair
    augment_parser(['annotations'], required[SUBCOMMAND.PAIR], optional[SUBCOMMAND.PAIR])
    augment_parser(['max_proximity'] + list(PAIRING_DEFAULTS.keys()), optional[SUBCOMMAND.PAIR])

    # summary
    augment_parser(
        ['annotations', 'flanking_call_distance', 'split_call_distance', 'contig_call_distance', 'spanning_call_distance'],
        required[SUBCOMMAND.SUMMARY]
    )
    augment_parser(SUMMARY_DEFAULTS.keys(), optional[SUBCOMMAND.SUMMARY])
    augment_parser(['dgv_annotation'], optional[SUBCOMMAND.SUMMARY])

    # overlay arguments
    required[SUBCOMMAND.OVERLAY].add_argument('gene_name', help='Gene ID or gene alias to be drawn')
    augment_parser(['annotations'], required[SUBCOMMAND.OVERLAY])
    augment_parser(['drawing_width_iter_increase', 'max_drawing_retries', 'width', 'min_mapping_quality'], optional[SUBCOMMAND.OVERLAY])
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--buffer_length', default=0, type=int, help='minimum genomic length to plot on either side of the target gene')
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--read_depth_plot', dest='read_depth_plots', metavar='<axis name STR> <bam FILEPATH> [density FLOAT] [ymax INT] [stranded BOOL]',
        nmin=2, nmax=5, help='bam file to use as data for plotting read_depth', action=RangeAppendAction)
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--marker', dest='markers', metavar='<label STR> <start INT> [end INT]', nmin=2, nmax=3,
        help='Marker on the diagram given by genomic position, May be a single position or a range. '
        'The label should be a short descriptor to avoid overlapping labels on the diagram',
        action=RangeAppendAction)

    args = MavisNamespace(**parser.parse_args().__dict__)
    if args.command == SUBCOMMAND.OVERLAY:
        args = check_overlay_args(args, parser)

    if args.command == SUBCOMMAND.VALIDATE:
        args.aligner_version = get_aligner_version(args.aligner)

    LOG('MAVIS: {}'.format(__version__))
    LOG('hostname:', platform.node(), time_stamp=False)
    log_arguments(args)
    rfile_args = args

    if args.command == SUBCOMMAND.PIPELINE:  # load the configuration file
        config = MavisConfig.read(args.config)
        config.output = args.output
        config.skip_stage = args.skip_stage
        config.command = SUBCOMMAND.PIPELINE
        rfile_args = config.reference
        args = config

    # try checking the input files exist
    try:
        inputs = []
        for fname in args.inputs:
            expanded = bash_expands(fname)
            inputs.extend(expanded)
            if not expanded:
                parser.error('--inputs file(s) {} do not exist'.format(args.inputs))
        args.inputs = set([os.path.abspath(f) for f in inputs])
    except AttributeError:
        pass

    # convert reference files to objects to store both content and name for rewrite
    for arg, loader in [
        ('reference_genome', _annotate.file_io.load_reference_genome),
        ('masking', _annotate.file_io.load_masking_regions),
        ('template_metadata', _annotate.file_io.load_templates),
        ('dgv_annotation', _annotate.file_io.load_masking_regions),
        ('annotations', _annotate.file_io.load_annotations)
    ]:
        if arg in rfile_args:
            rfile_args[arg] = _annotate.file_io.ReferenceFile(loader, *rfile_args[arg])

    # throw an error if MAVIS can't find the aligner reference
    if rfile_args.get('aligner_reference', False):
        if not rfile_args.aligner_reference:
            parser.error('--aligner_reference file does not exist at: {}'.format(rfile_args.aligner_reference))

    # for specific cases throw an argument error if missing annotations
    if any([
        args.command == SUBCOMMAND.CLUSTER and args.uninformative_filter,
        args.command == SUBCOMMAND.CONFIG and any([PROTOCOL.TRANS in values for values in args.library]) and SUBCOMMAND.VALIDATE not in args.skip_stage,
        args.command == SUBCOMMAND.VALIDATE and args.protocol == PROTOCOL.TRANS,
        args.command in {SUBCOMMAND.PAIR, SUBCOMMAND.ANNOTATE, SUBCOMMAND.SUMMARY, SUBCOMMAND.OVERLAY, SUBCOMMAND.PIPELINE},
    ]):
        try:
            rfile_args.annotations.files_exist()
        except FileNotFoundError:
            parser.error('--annotations file(s) are required and do not exist')

    # decide which main function to execute
    if args.command == SUBCOMMAND.CLUSTER:
        cluster_main.main(**args, start_time=start_time)
    elif args.command == SUBCOMMAND.VALIDATE:
        validate_main.main(**args, start_time=start_time)
    elif args.command == SUBCOMMAND.ANNOTATE:
        annotate_main.main(**args, start_time=start_time)
    elif args.command == SUBCOMMAND.PAIR:
        pairing_main.main(**args, start_time=start_time)
    elif args.command == SUBCOMMAND.SUMMARY:
        summary_main.main(**args, start_time=start_time)
    elif args.command == SUBCOMMAND.CONVERT:
        del args.command
        convert_main(**args)
    elif args.command == SUBCOMMAND.OVERLAY:
        del args.command
        overlay_main(**args)
    elif args.command == SUBCOMMAND.CONFIG:
        generate_config(args, parser, log=LOG)
    elif args.command == SUBCOMMAND.CHECKER:
        pass#return EXIT_OK if check_completion(args.output) else EXIT_ERROR
    elif args.command == SUBCOMMAND.SCHEDULE:
        build_file = os.path.join(args.output, 'build.cfg')
        pipeline = _schedule.pipeline.Pipeline.read_build_file(build_file)
        try:
            pipeline.check_status(log=LOG, submit=args.submit)
        finally:
            LOG('rewriting:', build_file)
            pipeline.write_build_file(build_file)
    else:  # PIPELINE
        config.reference = rfile_args
        pipeline = _schedule.pipeline.Pipeline.build(config)
        build_file = os.path.join(config.output, 'build.cfg')
        LOG('writing:', build_file)
        pipeline.write_build_file(build_file)

    duration = int(time.time()) - start_time
    hours = duration - duration % 3600
    minutes = duration - hours - (duration - hours) % 60
    seconds = duration - hours - minutes
    LOG(
        'run time (hh/mm/ss): {}:{:02d}:{:02d}'.format(hours // 3600, minutes // 60, seconds),
        time_stamp=False)
    LOG('run time (s): {}'.format(duration), time_stamp=False)
    return EXIT_OK


if __name__ == '__main__':
    exit(main())
