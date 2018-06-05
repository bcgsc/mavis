#!python
import argparse
import logging
import platform
import os
import time
import sys

import tab

from . import __version__
from .align import get_aligner_version
from . import annotate as _annotate
from .annotate import main as annotate_main
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .cluster import main as cluster_main
from . import config as _config
from .constants import SUBCOMMAND, PROTOCOL, float_fraction, EXIT_OK
from .error import DrawingFitError
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS, DiagramSettings
from .illustrate.diagram import draw_multi_transcript_overlay
from .illustrate.scatter import bam_to_scatter
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .pairing import main as pairing_main
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .summary import main as summary_main
from .tools import convert_tool_output, SUPPORTED_TOOL
from . import util as _util
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from .validate import main as validate_main
from .schedule import pipeline as _pipeline


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
                _util.LOG('Found target gene: {}(aka. {}) {}:{}-{}'.format(gene.name, gene.aliases, gene.chr, gene.start, gene.end))
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
                log=_util.LOG
            )
            break
        except DrawingFitError as err:
            if attempts > max_drawing_retries:
                raise err
            _util.LOG('Drawing fit: extending window', drawing_width_iter_increase)
            settings.width += drawing_width_iter_increase
            attempts += 1

    svg_output_file = os.path.join(output, '{}_{}_overlay.svg'.format(gene_to_draw.name, gene_name))
    _util.LOG('writing:', svg_output_file)

    canvas.saveas(svg_output_file)


def convert_main(inputs, outputfile, file_type, strand_specific=False, assume_no_untemplated=True):
    bpp_results = convert_tool_output(inputs, file_type, strand_specific, _util.LOG, True, assume_no_untemplated=assume_no_untemplated)
    if os.path.dirname(outputfile):
        _util.mkdirp(os.path.dirname(outputfile))
    _util.output_tabbed_file(bpp_results, outputfile)


def main(argv=None):
    """
    sets up the parser and checks the validity of command line args
    loads reference files and redirects into subcommand main functions

    Args:
        argv (list): List of arguments, defaults to command line arguments
    """
    if argv is None:  # need to do at run time or patching will not behave as expected
        argv = sys.argv[1:]
    start_time = int(time.time())

    parser = argparse.ArgumentParser(formatter_class=_config.CustomHelpFormatter)
    _config.augment_parser(['version'], parser)
    subp = parser.add_subparsers(dest='command', help='specifies which step/stage in the pipeline or which subprogram to use')
    subp.required = True
    required = {}  # hold required argument group by subparser command name
    optional = {}  # hold optional argument group by subparser command name
    for command in SUBCOMMAND.values():
        subparser = subp.add_parser(command, formatter_class=_config.CustomHelpFormatter, add_help=False)
        required[command] = subparser.add_argument_group('required arguments')
        optional[command] = subparser.add_argument_group('optional arguments')
        _config.augment_parser(['help', 'version', 'log', 'log_level'], optional[command])

    # config arguments
    required[SUBCOMMAND.CONFIG].add_argument('-w', '--write', help='path to the new configuration file', required=True, metavar='FILEPATH')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--library',
        metavar='<name> {genome,transcriptome} {diseased,normal} [strand_specific] [/path/to/bam/file]',
        action=_config.RangeAppendAction, help='configuration for libraries to be analyzed by mavis', nmin=3, nmax=5
    )
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--input', help='path to an input file or filter for mavis followed by the library names it '
        'should be used for', nmin=2, action=_config.RangeAppendAction, metavar='FILEPATH <name> [<name> ...]'
    )
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--assign', help='library name followed by path(s) to input file(s) or filter names. This represents the list'
        ' of inputs that should be used for the library', action=_config.RangeAppendAction, nmin=2,
        metavar='<name> FILEPATH [FILEPATH ...]')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--genome_bins', default=_util.get_env_variable('genome_bins', 100), type=int, metavar=_config.get_metavar(int),
        help='number of bins/samples to use in calculating the fragment size stats for genomes')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--transcriptome_bins', default=_util.get_env_variable('transcriptome_bins', 500), type=int, metavar=_config.get_metavar(int),
        help='number of genes to use in calculating the fragment size stats for genomes')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--distribution_fraction', default=_util.get_env_variable('distribution_fraction', 0.97), type=float_fraction, metavar=_config.get_metavar(float),
        help='the proportion of the distribution of calculated fragment sizes to use in determining the stdev')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--convert', nmin=3,
        metavar='<alias> FILEPATH [FILEPATH ...] {{{}}} [stranded]'.format(','.join(SUPPORTED_TOOL.values())),
        help='input file conversion for internally supported tools', action=_config.RangeAppendAction)
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--external_conversion', metavar=('<alias>', '<"command">'), nargs=2, default=[],
        help='alias for use in inputs and full command (quoted)', action='append')
    optional[SUBCOMMAND.CONFIG].add_argument(
        '--add_defaults', default=False, action='store_true', help='write current defaults for all non-specified options to the config output')
    _config.augment_parser(['annotations'], optional[SUBCOMMAND.CONFIG])
    # add the optional annotations file (only need this is auto generating bam stats for the transcriptome)
    _config.augment_parser(['skip_stage'], optional[SUBCOMMAND.CONFIG])

    # convert
    required[SUBCOMMAND.CONVERT].add_argument(
        '--file_type', choices=sorted(SUPPORTED_TOOL.values()),
        required=True, help='Indicates the input file type to be parsed')
    _config.augment_parser(['strand_specific', 'assume_no_untemplated'], optional[SUBCOMMAND.CONVERT])
    required[SUBCOMMAND.CONVERT].add_argument('--outputfile', '-o', required=True, help='path to the outputfile', metavar='FILEPATH')

    for command in set(SUBCOMMAND.values()) - {SUBCOMMAND.CONFIG, SUBCOMMAND.CONVERT}:
        required[command].add_argument('-o', '--output', help='path to the output directory', required=True)

    # pipeline
    _config.augment_parser(['config'], required[SUBCOMMAND.SETUP])
    optional[SUBCOMMAND.SETUP].add_argument(
        '--skip_stage', choices=[SUBCOMMAND.CLUSTER, SUBCOMMAND.VALIDATE], action='append', default=[],
        help='Use flag once per stage to skip. Can skip clustering or validation or both')

    # schedule arguments
    optional[SUBCOMMAND.SCHEDULE].add_argument('--submit', action='store_true', default=False, help='submit jobs to the the scheduler specified')
    optional[SUBCOMMAND.SCHEDULE].add_argument('--resubmit', action='store_true', default=False, help='resubmit jobs in error states to the the scheduler specified')

    # add the inputs argument
    for command in [SUBCOMMAND.CLUSTER, SUBCOMMAND.ANNOTATE, SUBCOMMAND.VALIDATE, SUBCOMMAND.PAIR, SUBCOMMAND.SUMMARY, SUBCOMMAND.CONVERT]:
        required[command].add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')

    # cluster
    _config.augment_parser(
        ['library', 'protocol', 'strand_specific', 'disease_status'],
        required[SUBCOMMAND.CLUSTER])
    _config.augment_parser(list(CLUSTER_DEFAULTS.keys()) + ['masking', 'annotations'], optional[SUBCOMMAND.CLUSTER])
    optional[SUBCOMMAND.CLUSTER].add_argument('--batch_id', help='batch id to use for prefix of split files', type=_config.nameable_string)
    optional[SUBCOMMAND.CLUSTER].add_argument('--split_only', help='Cluster the files or simply split them without clustering', type=tab.cast_boolean)

    # validate
    _config.augment_parser(
        [
            'library', 'protocol', 'bam_file', 'read_length', 'stdev_fragment_size',
            'median_fragment_size', 'strand_specific', 'reference_genome', 'aligner_reference'
        ],
        required[SUBCOMMAND.VALIDATE]
    )
    _config.augment_parser(VALIDATION_DEFAULTS.keys(), optional[SUBCOMMAND.VALIDATE])
    _config.augment_parser(['masking', 'annotations'], optional[SUBCOMMAND.VALIDATE])

    # annotate
    _config.augment_parser(
        ['library', 'protocol', 'annotations', 'reference_genome'],
        required[SUBCOMMAND.ANNOTATE]
    )
    _config.augment_parser(['max_proximity', 'masking', 'template_metadata'], optional[SUBCOMMAND.ANNOTATE])
    _config.augment_parser(list(_annotate.constants.DEFAULTS.keys()) + list(ILLUSTRATION_DEFAULTS.keys()), optional[SUBCOMMAND.ANNOTATE])

    # pair
    _config.augment_parser(['annotations'], required[SUBCOMMAND.PAIR], optional[SUBCOMMAND.PAIR])
    _config.augment_parser(['max_proximity'] + list(PAIRING_DEFAULTS.keys()), optional[SUBCOMMAND.PAIR])

    # summary
    _config.augment_parser(
        ['annotations', 'flanking_call_distance', 'split_call_distance', 'contig_call_distance', 'spanning_call_distance'],
        required[SUBCOMMAND.SUMMARY]
    )
    _config.augment_parser(SUMMARY_DEFAULTS.keys(), optional[SUBCOMMAND.SUMMARY])
    _config.augment_parser(['dgv_annotation'], optional[SUBCOMMAND.SUMMARY])

    # overlay arguments
    required[SUBCOMMAND.OVERLAY].add_argument('gene_name', help='Gene ID or gene alias to be drawn')
    _config.augment_parser(['annotations'], required[SUBCOMMAND.OVERLAY])
    _config.augment_parser(['drawing_width_iter_increase', 'max_drawing_retries', 'width', 'min_mapping_quality'], optional[SUBCOMMAND.OVERLAY])
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--buffer_length', default=0, type=int, help='minimum genomic length to plot on either side of the target gene')
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--read_depth_plot', dest='read_depth_plots', metavar='<axis name STR> <bam FILEPATH> [density FLOAT] [ymax INT] [stranded BOOL]',
        nmin=2, nmax=5, help='bam file to use as data for plotting read_depth', action=_config.RangeAppendAction)
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--marker', dest='markers', metavar='<label STR> <start INT> [end INT]', nmin=2, nmax=3,
        help='Marker on the diagram given by genomic position, May be a single position or a range. '
        'The label should be a short descriptor to avoid overlapping labels on the diagram',
        action=_config.RangeAppendAction)

    args = _util.MavisNamespace(**parser.parse_args(argv).__dict__)
    if args.command == SUBCOMMAND.OVERLAY:
        args = check_overlay_args(args, parser)

    if args.command == SUBCOMMAND.VALIDATE:
        args.aligner_version = get_aligner_version(args.aligner)

    log_conf = {'format': '{message}', 'style': '{', 'level': args.log_level}

    original_logging_handlers = logging.root.handlers[:]
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)
    if args.log:  # redirect stdout AND stderr to a log file
        log_conf['filename'] = args.log
    logging.basicConfig(**log_conf)

    _util.LOG('MAVIS: {}'.format(__version__))
    _util.LOG('hostname:', platform.node(), time_stamp=False)
    _util.log_arguments(args)
    rfile_args = args

    if args.command == SUBCOMMAND.SETUP:  # load the configuration file
        config = _config.MavisConfig.read(args.config)
        config.output = args.output
        config.skip_stage = args.skip_stage
        config.command = SUBCOMMAND.SETUP
        rfile_args = config.reference
        args = config

    # try checking the input files exist
    try:
        args.inputs = _util.bash_expands(*args.inputs)
    except AttributeError:
        pass
    except FileNotFoundError:
        parser.error('--inputs file(s) for {} {} do not exist'.format(args.command, args.inputs))

    # convert reference files to objects to store both content and name for rewrite
    for arg in [f for f in _annotate.file_io.REFERENCE_DEFAULTS.keys() if f != 'aligner_reference']:
        try:
            rfile_args[arg] = _annotate.file_io.ReferenceFile(arg, assert_exists=True, *rfile_args[arg])
        except AttributeError:
            pass
        except FileNotFoundError:
            parser.error('--{} The file specified does not exist: {}'.format(arg, rfile_args[arg]))

    # throw an error if MAVIS can't find the aligner reference
    try:
        rfile_args.aligner_reference = _annotate.file_io.ReferenceFile('aligner_reference', rfile_args.aligner_reference, assert_exists=True)
    except AttributeError:
        pass
    except FileNotFoundError:
        parser.error('--aligner_reference file does not exist at: {}'.format(rfile_args.aligner_reference))

    # for specific cases throw an argument error if missing annotations
    if any([
        args.command == SUBCOMMAND.CLUSTER and args.uninformative_filter,
        args.command == SUBCOMMAND.CONFIG and any([PROTOCOL.TRANS in values for values in args.library]) and SUBCOMMAND.VALIDATE not in args.skip_stage,
        args.command == SUBCOMMAND.VALIDATE and args.protocol == PROTOCOL.TRANS,
        args.command in {SUBCOMMAND.PAIR, SUBCOMMAND.ANNOTATE, SUBCOMMAND.SUMMARY, SUBCOMMAND.OVERLAY, SUBCOMMAND.SETUP},
    ]):
        try:
            rfile_args.annotations.files_exist(not_empty=True)
        except FileNotFoundError:
            parser.error('--annotations file(s) are required and do not exist')

    # decide which main function to execute
    ret_val = EXIT_OK
    command = args.command
    log_to_file = args.get('log', None)

    # discard any arguments needed for redirect/setup only
    for init_arg in ['command', 'log', 'log_level']:
        args.discard(init_arg)

    try:
        if command == SUBCOMMAND.CLUSTER:
            ret_val = cluster_main.main(**args, start_time=start_time)
        elif command == SUBCOMMAND.VALIDATE:
            validate_main.main(**args, start_time=start_time)
        elif command == SUBCOMMAND.ANNOTATE:
            annotate_main.main(**args, start_time=start_time)
        elif command == SUBCOMMAND.PAIR:
            pairing_main.main(**args, start_time=start_time)
        elif command == SUBCOMMAND.SUMMARY:
            summary_main.main(**args, start_time=start_time)
        elif command == SUBCOMMAND.CONVERT:
            convert_main(**args)
        elif command == SUBCOMMAND.OVERLAY:
            overlay_main(**args)
        elif command == SUBCOMMAND.CONFIG:
            _config.generate_config(args, parser, log=_util.LOG)
        elif command == SUBCOMMAND.SCHEDULE:
            build_file = os.path.join(args.output, 'build.cfg')
            args.discard('output')
            pipeline = _pipeline.Pipeline.read_build_file(build_file)
            try:
                code = pipeline.check_status(log=_util.LOG, **args)
            finally:
                _util.LOG('rewriting:', build_file)
                pipeline.write_build_file(build_file)
            if code != EXIT_OK:
                sys.exit(code)  # EXIT
        else:  # PIPELINE
            config.reference = rfile_args
            pipeline = _pipeline.Pipeline.build(config)
            build_file = os.path.join(config.output, 'build.cfg')
            _util.LOG('writing:', build_file)
            pipeline.write_build_file(build_file)

        duration = int(time.time()) - start_time
        hours = duration - duration % 3600
        minutes = duration - hours - (duration - hours) % 60
        seconds = duration - hours - minutes
        _util.LOG(
            'run time (hh/mm/ss): {}:{:02d}:{:02d}'.format(hours // 3600, minutes // 60, seconds),
            time_stamp=False)
        _util.LOG('run time (s): {}'.format(duration), time_stamp=False)
        return ret_val
    except Exception as err:
        if log_to_file:
            logging.exception(err)  # capture the error in the logging output file
        raise err
    finally:
        for handler in logging.root.handlers:
            logging.root.removeHandler(handler)
        for handler in original_logging_handlers:
            logging.root.addHandler(handler)


if __name__ == '__main__':
    main()
