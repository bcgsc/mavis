#!python
import argparse
import json
import logging
import os
import platform
import sys
import time
from typing import Dict, List, Optional

from mavis_config import validate_config
from mavis_config.constants import SUBCOMMAND

from . import __version__
from . import config as _config
from . import util as _util
from .align import get_aligner_version
from .annotate import main as annotate_main
from .cluster import main as cluster_main
from .convert import SUPPORTED_TOOL, convert_tool_output
from .overlay import check_overlay_args
from .overlay import main as overlay_main
from .pairing import main as pairing_main
from .summary import main as summary_main
from .util import filepath
from .validate import main as validate_main


def convert_main(inputs, outputfile, file_type, strand_specific=False, assume_no_untemplated=True):
    bpp_results = convert_tool_output(
        inputs,
        file_type,
        strand_specific,
        True,
        assume_no_untemplated=assume_no_untemplated,
    )
    if os.path.dirname(outputfile):
        _util.mkdirp(os.path.dirname(outputfile))
    _util.output_tabbed_file(bpp_results, outputfile)


def create_parser(argv):
    parser = argparse.ArgumentParser(formatter_class=_config.CustomHelpFormatter)
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s version ' + __version__,
        help='Outputs the version number',
    )
    subp = parser.add_subparsers(
        dest='command', help='specifies which step/stage in the pipeline or which subprogram to use'
    )
    subp.required = True
    required = {}  # hold required argument group by subparser command name
    optional = {}  # hold optional argument group by subparser command name
    for command in SUBCOMMAND.values():
        subparser = subp.add_parser(
            command, formatter_class=_config.CustomHelpFormatter, add_help=False
        )
        required[command] = subparser.add_argument_group('required arguments')
        optional[command] = subparser.add_argument_group('optional arguments')
        optional[command].add_argument(
            '-h', '--help', action='help', help='show this help message and exit'
        )
        optional[command].add_argument(
            '-v',
            '--version',
            action='version',
            version='%(prog)s version ' + __version__,
            help='Outputs the version number',
        )
        optional[command].add_argument('--log', help='redirect stdout to a log file', default=None)
        optional[command].add_argument(
            '--log_level',
            help='level of logging to output',
            choices=['INFO', 'DEBUG'],
            default='INFO',
        )
        if command not in SUBCOMMAND.CONVERT:
            optional[command].add_argument(
                '--config', '-c', help='path to the JSON config file', type=filepath, required=True
            )

    # convert
    required[SUBCOMMAND.CONVERT].add_argument(
        '--file_type',
        choices=sorted(SUPPORTED_TOOL.values()),
        required=True,
        help='Indicates the input file type to be parsed',
    )
    optional[SUBCOMMAND.CONVERT].add_argument(
        '--strand_specific', type=_util.cast_boolean, default=False
    )
    optional[SUBCOMMAND.CONVERT].add_argument(
        '--assume_no_untemplated', type=_util.cast_boolean, default=True
    )
    for command in [SUBCOMMAND.CONVERT, SUBCOMMAND.SETUP]:
        required[command].add_argument(
            '--outputfile', '-o', required=True, help='path to the outputfile', metavar='FILEPATH'
        )

    for command in set(SUBCOMMAND.values()) - {SUBCOMMAND.CONVERT, SUBCOMMAND.SETUP}:
        required[command].add_argument(
            '-o', '--output', help='path to the output directory', required=True
        )

    # add the inputs argument
    for command in [
        SUBCOMMAND.CLUSTER,
        SUBCOMMAND.ANNOTATE,
        SUBCOMMAND.VALIDATE,
        SUBCOMMAND.PAIR,
        SUBCOMMAND.SUMMARY,
        SUBCOMMAND.CONVERT,
    ]:
        required[command].add_argument(
            '-n',
            '--inputs',
            nargs='+',
            help='path to the input files',
            required=True,
            metavar='FILEPATH',
        )

    # library specific commands
    for command in [SUBCOMMAND.CLUSTER, SUBCOMMAND.VALIDATE, SUBCOMMAND.ANNOTATE]:
        required[command].add_argument(
            '--library', '-l', required=True, help='The library to run the current step on'
        )

    # overlay arguments
    required[SUBCOMMAND.OVERLAY].add_argument('gene_name', help='Gene ID or gene alias to be drawn')
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--buffer_length',
        default=0,
        type=int,
        help='minimum genomic length to plot on either side of the target gene',
    )
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--read_depth_plot',
        dest='read_depth_plots',
        metavar='<axis name STR> <bam FILEPATH> [density FLOAT] [ymax INT] [stranded BOOL]',
        nmin=2,
        nmax=5,
        help='bam file to use as data for plotting read_depth',
        action=_config.RangeAppendAction,
    )
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--marker',
        dest='markers',
        metavar='<label STR> <start INT> [end INT]',
        nmin=2,
        nmax=3,
        help='Marker on the diagram given by genomic position, May be a single position or a range. '
        'The label should be a short descriptor to avoid overlapping labels on the diagram',
        action=_config.RangeAppendAction,
    )

    return parser, parser.parse_args(argv)


def main(argv: Optional[List[str]] = None):
    """
    sets up the parser and checks the validity of command line args
    loads reference files and redirects into subcommand main functions

    Args:
        argv: List of arguments, defaults to command line arguments
    """
    if argv is None:  # need to do at run time or patching will not behave as expected
        argv = sys.argv[1:]
    start_time = int(time.time())
    parser, args = create_parser(argv)

    if args.command == SUBCOMMAND.OVERLAY:
        args = check_overlay_args(args, parser)

    log_conf = {
        'format': '{asctime} [{levelname}] {message}',
        'style': '{',
        'level': args.log_level,
    }

    original_logging_handlers = logging.root.handlers[:]
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)
    if args.log:  # redirect stdout AND stderr to a log file
        log_conf['filename'] = args.log
    logging.basicConfig(**log_conf)

    _util.logger.info(f'MAVIS: {__version__}')
    _util.logger.info(f'hostname: {platform.node()}')
    _util.log_arguments(args)

    config: Dict = dict()

    try:
        if args.command != SUBCOMMAND.CONVERT:
            with open(args.config, 'r') as fh:
                config = json.load(fh)
                validate_config(
                    config,
                    args.command,
                )
    except AttributeError as err:
        raise err

    if args.command == SUBCOMMAND.VALIDATE:
        args.aligner_version = get_aligner_version(config['validate.aligner'])
    # try checking the input files exist
    try:
        args.inputs = _util.bash_expands(*args.inputs)
    except AttributeError:
        pass
    except FileNotFoundError:
        parser.error('--inputs file(s) for {} {} do not exist'.format(args.command, args.inputs))

    # decide which main function to execute
    command = args.command

    try:
        if command == SUBCOMMAND.CLUSTER:
            cluster_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
                library=args.library,
            )
        elif command == SUBCOMMAND.VALIDATE:
            validate_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
                library=args.library,
            )
        elif command == SUBCOMMAND.ANNOTATE:
            annotate_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
                library=args.library,
            )
        elif command == SUBCOMMAND.PAIR:
            pairing_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
            )
        elif command == SUBCOMMAND.SUMMARY:
            summary_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
            )
        elif command == SUBCOMMAND.CONVERT:
            convert_main(
                args.inputs,
                args.outputfile,
                args.file_type,
                args.strand_specific,
                args.assume_no_untemplated,
            )
        elif command == SUBCOMMAND.SETUP:
            # add bam stats to the config if missing
            if not config.get('skip_stage.validate'):
                _config.add_bamstats_to_config(config)
            _util.logger.info(f'writing: {args.outputfile}')
            with open(args.outputfile, 'w') as fh:
                fh.write(json.dumps(config, sort_keys=True, indent='  '))
        else:
            overlay_main(
                buffer_length=args.buffer_length,
                gene_name=args.gene_name,
                markers=args.markers,
                read_depth_plots=args.read_depth_plots,
                config=config,
                output=args.output,
            )

        duration = int(time.time()) - start_time
        hours = duration - duration % 3600
        minutes = duration - hours - (duration - hours) % 60
        seconds = duration - hours - minutes
        _util.logger.info(
            'run time (hh/mm/ss): {}:{:02d}:{:02d}'.format(hours // 3600, minutes // 60, seconds)
        )
        _util.logger.info(f'run time (s): {duration}')
    except Exception as err:
        raise err
    finally:
        try:
            for handler in logging.root.handlers:
                logging.root.removeHandler(handler)
            for handler in original_logging_handlers:
                logging.root.addHandler(handler)
        except Exception as err:
            print(err)


if __name__ == '__main__':
    main()
