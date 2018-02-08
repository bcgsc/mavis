#!python
import argparse
import platform
import os
import re
import subprocess
import sys
import time

import tab

from . import __version__
from .annotate.base import BioInterval
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .annotate.file_io import load_annotations, load_masking_regions, load_reference_genome, load_templates
from .annotate.main import main as annotate_main
from .bam.read import get_samtools_version
from .blat import get_blat_version
from .checker import check_completion
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .cluster.main import main as cluster_main
from .config import augment_parser, MavisConfig, generate_config, CustomHelpFormatter, RangeAppendAction
from .constants import SUBCOMMAND, PROTOCOL
from .error import DrawingFitError
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS, DiagramSettings
from .illustrate.diagram import draw_multi_transcript_overlay
from .illustrate.scatter import bam_to_scatter
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .pairing.main import main as pairing_main
from .submit import SubmissionScript, SCHEDULER_CONFIG
from .submit import STD_OPTIONS as STD_SUBMIT_OPTIONS
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .summary.main import main as summary_main
from .tools import convert_tool_output, SUPPORTED_TOOL
from .util import bash_expands, log, log_arguments, MavisNamespace, mkdirp, output_tabbed_file
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from .validate.main import main as validate_main


VALIDATION_PASS_PATTERN = '*.validation-passed.tab'
ANNOTATION_PASS_PATTERN = 'annotations.tab'
PROGNAME = 'mavis'
EXIT_OK = 0
EXIT_ERROR = 1


def build_validate_command(config, libconf, inputfile, outputdir):
    """
    returns the mavis command for running the validate step
    """
    args = {
        'masking': config.reference.masking_filename,
        'reference_genome': config.reference.reference_genome_filename,
        'aligner_reference': config.reference.aligner_reference,
        'annotations': config.reference.annotations_filename,
        'library': libconf.library,
        'bam_file': libconf.bam_file,
        'protocol': libconf.protocol,
        'read_length': libconf.read_length,
        'stdev_fragment_size': libconf.stdev_fragment_size,
        'median_fragment_size': libconf.median_fragment_size,
        'strand_specific': libconf.strand_specific
    }
    args.update(config.validate.items())
    args.update({k: v for k, v in libconf.items() if k in args})

    command = ['{} {}'.format(PROGNAME, SUBCOMMAND.VALIDATE)]
    for argname, value in args.items():
        if isinstance(value, str):
            command.append('--{} "{}"'.format(argname, value))
        else:
            command.append('--{} {}'.format(argname, value))
    command.append('--input {}'.format(inputfile))
    command.append('--output {}'.format(outputdir))
    return ' \\\n\t'.join(command) + '\n'


def build_annotate_command(config, libconf, inputfile, outputdir):
    """
    returns the mavis command for running the annotate step
    """
    args = {
        'reference_genome': config.reference.reference_genome_filename,
        'annotations': config.reference.annotations_filename,
        'template_metadata': config.reference.template_metadata_filename,
        'masking': config.reference.masking_filename,
        'min_orf_size': config.annotate.min_orf_size,
        'max_orf_cap': config.annotate.max_orf_cap,
        'library': libconf.library,
        'protocol': libconf.protocol,
        'min_domain_mapping_match': config.annotate.min_domain_mapping_match,
        'domain_name_regex_filter': config.illustrate.domain_name_regex_filter,
        'max_proximity': config.cluster.max_proximity
    }
    args.update(config.annotate.items())
    args.update({k: v for k, v in libconf.items() if k in args})
    command = ['{} {}'.format(PROGNAME, SUBCOMMAND.ANNOTATE)]
    for argname, value in args.items():
        if isinstance(value, str):
            command.append('--{} "{}"'.format(argname, value))
        else:
            command.append('--{} {}'.format(argname, value))
    command.append('--inputs {}'.format(inputfile))
    command.append('--output {}'.format(outputdir))
    return ' \\\n\t'.join(command) + '\n'


def run_conversion(config, libconf, conversion_dir, assume_no_untemplated=True):
    """
    Converts files if not already converted. Returns a list of filenames
    """
    inputs = []
    # run the conversions
    for input_file in libconf.inputs:
        output_filename = os.path.join(conversion_dir, input_file + '.tab')
        if input_file in config.convert:
            if not os.path.exists(output_filename):
                command = config.convert[input_file]
                if command[0] == 'convert_tool_output':
                    log('converting input command:', command)
                    output_tabbed_file(convert_tool_output(*command[1:], log=log, assume_no_untemplated=assume_no_untemplated), output_filename)
                else:
                    command = ' '.join(command) + ' -o {}'.format(output_filename)
                    log('converting input command:')
                    log('>>>', command, time_stamp=False)
                    subprocess.check_output(command, shell=True)
            inputs.append(output_filename)
        else:
            inputs.append(input_file)
    return inputs


def main_pipeline(config):
    """
    runs the pipeline subcommand. Runs clustering (or just splitting if clustering is skipped) and sets up
    submission scripts for the other pipeline stages/steps
    """
    from shortuuid import uuid
    batch_id = 'batch-' + str(uuid())
    conversion_dir = mkdirp(os.path.join(config.output, 'converted_inputs'))
    pairing_inputs = []
    submitall = []
    jobid_var_index = 0
    config.output = os.path.abspath(config.output)

    def get_prefix(filename):
        return re.sub(r'\.[^\.]+$', '', os.path.basename(filename))

    job_name_by_output = {}

    for libconf in config.libraries.values():
        base = os.path.join(config.output, '{}_{}_{}'.format(libconf.library, libconf.disease_status, libconf.protocol))
        log('setting up the directory structure for', libconf.library, 'as', base)
        libconf.inputs = run_conversion(config, libconf, conversion_dir)

        # run the cluster stage
        cluster_output = mkdirp(os.path.join(base, SUBCOMMAND.CLUSTER))  # creates the output dir
        merge_args = {'batch_id': batch_id, 'output': cluster_output}
        merge_args['split_only'] = SUBCOMMAND.CLUSTER in config.skip_stage
        merge_args.update(config.reference.items())
        merge_args.update(config.cluster.items())
        merge_args.update(libconf.items())
        log('clustering', '(split only)' if merge_args['split_only'] else '')
        inputs = cluster_main(log_args=True, **merge_args)

        for inputfile in inputs:
            prefix = get_prefix(inputfile)  # will be batch id + job number
            dependency = ''
            if SUBCOMMAND.VALIDATE not in config.skip_stage:
                outputdir = mkdirp(os.path.join(base, SUBCOMMAND.VALIDATE, prefix))
                command = build_validate_command(config, libconf, inputfile, outputdir)
                # build the submission script
                options = {k: config.schedule[k] for k in STD_SUBMIT_OPTIONS}
                options['stdout'] = outputdir
                options['jobname'] = 'MV_{}_{}'.format(libconf.library, prefix)

                if libconf.is_trans():
                    options['memory_limit'] = config.schedule.trans_validation_memory
                else:
                    options['memory_limit'] = config.schedule.validation_memory
                script = SubmissionScript(command, config.schedule.scheduler, **options)
                scriptname = script.write(os.path.join(outputdir, 'submit.sh'))

                submitall.append('vjob{}=$({} {})'.format(jobid_var_index, SCHEDULER_CONFIG[config.schedule.scheduler].submit, scriptname))
                # for setting up subsequent jobs and holds
                outputfile = os.path.join(outputdir, VALIDATION_PASS_PATTERN)
                job_name_by_output[outputfile] = options['jobname']
                inputfile = outputfile
                dependency = SCHEDULER_CONFIG[config.schedule.scheduler].dependency('${{vjob{}##* }}'.format(jobid_var_index))
            # annotation cannot be skipped
            outputdir = mkdirp(os.path.join(base, SUBCOMMAND.ANNOTATE, prefix))
            command = build_annotate_command(config, libconf, inputfile, outputdir)

            options = {k: config.schedule[k] for k in STD_SUBMIT_OPTIONS}
            options['stdout'] = outputdir
            options['jobname'] = 'MA_{}_{}'.format(libconf.library, prefix)
            options['memory_limit'] = config.schedule.annotation_memory
            script = SubmissionScript(command, config.schedule.scheduler, **options)
            scriptname = script.write(os.path.join(outputdir, 'submit.sh'))
            submitall.append('ajob{}=$({} {} {})'.format(
                jobid_var_index, SCHEDULER_CONFIG[config.schedule.scheduler].submit, dependency, scriptname))
            outputfile = os.path.join(outputdir, ANNOTATION_PASS_PATTERN)
            pairing_inputs.append(outputfile)
            job_name_by_output[outputfile] = options['jobname']
            jobid_var_index += 1

    # set up scripts for the pairing held on all of the annotation jobs
    outputdir = mkdirp(os.path.join(config.output, SUBCOMMAND.PAIR))
    args = config.pairing.flatten()
    args.update({
        'output': outputdir,
        'annotations': config.reference.annotations_filename
    })
    command = ['{} {}'.format(PROGNAME, SUBCOMMAND.PAIR)]
    for arg, value in sorted(args.items()):
        if isinstance(value, str):
            command.append('--{} "{}"'.format(arg, value))
        else:
            command.append('--{} {}'.format(arg, value))
    command.append('--inputs {}'.format(' \\\n\t'.join(pairing_inputs)))
    command = ' \\\n\t'.join(command)

    options = {k: config.schedule[k] for k in STD_SUBMIT_OPTIONS}
    options['stdout'] = outputdir
    options['jobname'] = 'MP_{}'.format(batch_id)
    script = SubmissionScript(command, config.schedule.scheduler, **options)
    scriptname = script.write(os.path.join(outputdir, 'submit.sh'))

    submitall.append('jobid=$({} {} {})'.format(
        SCHEDULER_CONFIG[config.schedule.scheduler].submit,
        SCHEDULER_CONFIG[config.schedule.scheduler].dependency(
            ':'.join(['${{ajob{}##* }}'.format(i) for i in range(0, jobid_var_index)])),
        scriptname))

    # set up scripts for the summary held on the pairing job
    outputdir = mkdirp(os.path.join(config.output, SUBCOMMAND.SUMMARY))
    args = dict(
        output=outputdir,
        flanking_call_distance=config.pairing.flanking_call_distance,
        split_call_distance=config.pairing.split_call_distance,
        contig_call_distance=config.pairing.contig_call_distance,
        spanning_call_distance=config.pairing.spanning_call_distance,
        dgv_annotation=config.reference.dgv_annotation_filename,
        annotations=config.reference.annotations_filename,
        inputs=os.path.join(config.output, 'pairing/mavis_paired*.tab')
    )
    args.update(config.summary.items())
    command = ['{} {}'.format(PROGNAME, SUBCOMMAND.SUMMARY)]
    for arg, value in sorted(args.items()):
        if isinstance(value, str):
            command.append('--{} "{}"'.format(arg, value))
        else:
            command.append('--{} {}'.format(arg, value))
    command = ' \\\n\t'.join(command)

    options = {k: config.schedule[k] for k in STD_SUBMIT_OPTIONS}
    options['stdout'] = outputdir
    options['jobname'] = 'MS_{}'.format(batch_id)
    script = SubmissionScript(command, config.schedule.scheduler, **options)
    scriptname = script.write(os.path.join(outputdir, 'submit.sh'))

    submitall.append('{} {} {}'.format(
        SCHEDULER_CONFIG[config.schedule.scheduler].submit,
        SCHEDULER_CONFIG[config.schedule.scheduler].dependency('${jobid##* }'),
        scriptname))

    # now write a script at the top level to submit all
    submitallfile = os.path.join(config.output, 'submit_pipeline_{}.sh'.format(batch_id))
    log('writing:', submitallfile)
    with open(submitallfile, 'w') as fh:
        for line in submitall:
            fh.write(line + '\n')


def parse_overlay_args(parser, required, optional):
    """
    parse the overlay options and check the formatting
    """
    required.add_argument('gene_name', help='Gene ID or gene alias to be drawn')
    augment_parser(['annotations'], required, optional)
    augment_parser(['drawing_width_iter_increase', 'max_drawing_retries', 'width', 'min_mapping_quality'], optional)
    optional.add_argument(
        '--buffer_length', default=0, type=int, help='minimum genomic length to plot on either side of the target gene')
    optional.add_argument(
        '--read_depth_plot', dest='read_depth_plots', metavar='<axis name STR> <bam FILEPATH> [bin size INT] [ymax INT] [stranded BOOL]',
        nmin=2, nmax=5, help='bam file to use as data for plotting read_depth', action=RangeAppendAction)
    optional.add_argument(
        '--marker', dest='markers', metavar='<label STR> <start INT> [end INT]', nmin=2, nmax=3,
        help='Marker on the diagram given by genomic position, May be a single position or a range. '
        'The label should be a short descriptor to avoid overlapping labels on the diagram',
        action=RangeAppendAction)
    args = MavisNamespace(**parser.parse_args().__dict__)

    # check complex options
    for marker in args.markers:
        if len(marker) < 3:
            marker.append(marker[-1])
        try:
            marker[1] = int(marker[1])
            marker[2] = int(marker[2])
        except ValueError:
            parser.error('argument --marker: start and end must be integers: {}'.format(marker))

    defaults = [None, None, 1, None, True]
    bam_file, bin_size, ymax, stranded = range(1, 5)

    for plot in args.read_depth_plots:
        for i, d in enumerate(defaults):
            if i >= len(plot):
                plot.append(d)
        if not os.path.exists(plot[bam_file]):
            parser.error('argument --read_depth_plots: the bam file given does not exist: {}'.format(plot[bam_file]))
        try:
            plot[bin_size] = int(plot[bin_size])
        except ValueError:
            parser.error('argument --read_depth_plots: bin size must be an integer: {}'.format(plot[bin_size]))
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
    annotations, annotations_filename,
    drawing_width_iter_increase, max_drawing_retries, min_mapping_quality,
    **kwargs
):
    """
    generates an overlay diagram
    """
    # check options formatting
    gene_to_draw = None

    for chrom in annotations:
        for gene in annotations[chrom]:
            if gene_name in gene.aliases or gene_name == gene.name:
                gene_to_draw = gene
                log('Found target gene: {}(aka. {}) {}:{}-{}'.format(gene.name, gene.aliases, gene.chr, gene.start, gene.end))
                break
    if gene_to_draw is None:
        raise KeyError('Could not find gene alias or id in annotations file', gene_name)

    settings = DiagramSettings(**kwargs)

    x_start = max(gene_to_draw.start - buffer_length, 1)
    x_end = gene_to_draw.end + buffer_length

    plots = []
    for axis_name, bam_file, bin_size, ymax, stranded in read_depth_plots:
        # one plot per bam
        plots.append(bam_to_scatter(
            bam_file, gene_to_draw.chr, x_start, x_end,
            strand=gene.get_strand() if stranded else None,
            ymax=ymax,
            bin_size=bin_size,
            axis_name=axis_name,
            min_mapping_quality=min_mapping_quality
        ))
        log('reading:', bam_file)

    for i, (marker_name, marker_start, marker_end) in enumerate(markers):
        markers[i] = BioInterval(gene_to_draw.chr, marker_start, marker_end, name=marker_name)

    canvas = None
    attempts = 1
    while True:
        try:
            canvas = draw_multi_transcript_overlay(settings, gene_to_draw, vmarkers=markers, plots=plots)
            break
        except DrawingFitError as err:
            if attempts > max_drawing_retries:
                raise err
            log('Drawing fit: extending window', drawing_width_iter_increase)
            settings.width += drawing_width_iter_increase
            attempts += 1

    svg_output_file = os.path.join(output, '{}_{}_overlay.svg'.format(gene_to_draw.name, gene_name))
    log('writing:', svg_output_file)

    canvas.saveas(svg_output_file)


def convert_main(inputs, outputfile, file_type, strand_specific=False, assume_no_untemplated=True):
    bpp_results = []
    for filename in inputs:
        bpp_results.extend(convert_tool_output(filename, file_type, strand_specific, log, True, assume_no_untemplated=assume_no_untemplated))
    if os.path.dirname(outputfile):
        mkdirp(os.path.dirname(outputfile))
    output_tabbed_file(bpp_results, outputfile)


def main():
    def usage(err=None, detail=False):
        umsg = '\nusage: {} {{{}}} [-h] [-v]'.format(PROGNAME, ','.join(sorted(SUBCOMMAND.values())))
        helpmenu = """
required arguments:

    pipeline_step
        specifies which step in the pipeline or which subprogram
        should be run. See possible input values above

optional arguments:
    -h, --help
        bring up this help menu
    -v, --version
        output the version number

To bring up individual help menus for a given pipeline step
use the -h/--help option

    >>> {} <pipeline step> -h
    """.format(PROGNAME)
        print(umsg)
        if detail:
            print(helpmenu)
        if err:
            print('{}: error:'.format(PROGNAME), err, '\n')
            return EXIT_ERROR
        return EXIT_OK

    start_time = int(time.time())

    if len(sys.argv) < 2:
        return usage('the <pipeline step> argument is required')
    elif sys.argv[1] in ['-h', '--help']:
        return usage(detail=True)
    elif sys.argv[1] in ['-v', '--version']:
        print('{} version {}'.format('mavis', __version__))
        return EXIT_OK

    pstep = sys.argv.pop(1)
    sys.argv[0] = '{} {}'.format(sys.argv[0], pstep)

    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter, add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    augment_parser(['help', 'version'], optional)

    if pstep == SUBCOMMAND.CONFIG:
        generate_config(parser, required, optional, log=log)
        return EXIT_OK
    elif pstep == SUBCOMMAND.CONVERT:
        required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
        required.add_argument(
            '--file_type', choices=sorted([t for t in SUPPORTED_TOOL.values() if t != 'mavis']),
            required=True, help='Indicates the input file type to be parsed')
        augment_parser(['strand_specific', 'assume_no_untemplated'], optional)
        required.add_argument('--outputfile', '-o', required=True, help='path to the outputfile', metavar='FILEPATH')
    else:
        required.add_argument('-o', '--output', help='path to the output directory', required=True)

        if pstep == SUBCOMMAND.PIPELINE:
            required.add_argument('config', help='path to the input pipeline configuration file', metavar='FILEPATH')
            optional.add_argument(
                '--skip_stage', choices=[SUBCOMMAND.CLUSTER, SUBCOMMAND.VALIDATE], action='append', default=[],
                help='Use flag once per stage to skip. Can skip clustering or validation or both')

        elif pstep == SUBCOMMAND.CLUSTER:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(['library', 'protocol', 'strand_specific', 'disease_status', 'annotations', 'masking'], required, optional)
            augment_parser(CLUSTER_DEFAULTS.keys(), optional)

        elif pstep == SUBCOMMAND.VALIDATE:
            required.add_argument('-n', '--input', help='path to the input file', required=True, metavar='FILEPATH')
            augment_parser(
                ['library', 'protocol', 'bam_file', 'read_length', 'stdev_fragment_size', 'median_fragment_size'] +
                ['strand_specific', 'annotations', 'reference_genome', 'aligner_reference', 'masking'],
                required, optional
            )
            augment_parser(VALIDATION_DEFAULTS.keys(), optional)

        elif pstep == SUBCOMMAND.ANNOTATE:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(
                ['library', 'protocol', 'annotations', 'reference_genome', 'masking', 'template_metadata'],
                required, optional
            )
            augment_parser(['max_proximity'], optional)
            augment_parser(list(ANNOTATION_DEFAULTS.keys()) + list(ILLUSTRATION_DEFAULTS.keys()), optional)

        elif pstep == SUBCOMMAND.PAIR:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(['annotations'], required, optional)
            augment_parser(['max_proximity'] + list(PAIRING_DEFAULTS.keys()), optional)

        elif pstep == SUBCOMMAND.SUMMARY:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(
                ['annotations', 'dgv_annotation', 'flanking_call_distance', 'split_call_distance', 'contig_call_distance', 'spanning_call_distance'],
                required, optional
            )
            augment_parser(SUMMARY_DEFAULTS.keys(), optional)

        elif pstep == SUBCOMMAND.CHECKER:

            args = parser.parse_args()
            success_flag = check_completion(args.output)
            return EXIT_OK if success_flag else EXIT_ERROR

        elif pstep != SUBCOMMAND.OVERLAY:

            raise NotImplementedError('invalid value for <pipeline step>', pstep)
    if pstep == SUBCOMMAND.OVERLAY:
        args = parse_overlay_args(parser, required, optional)
    else:
        args = MavisNamespace(**parser.parse_args().__dict__)

    if pstep == SUBCOMMAND.VALIDATE:
        args.samtools_version = get_samtools_version()
        args.blat_version = get_blat_version()

    log('MAVIS: {}'.format(__version__))
    log('hostname:', platform.node(), time_stamp=False)
    log_arguments(args)
    rargs = args

    if pstep == SUBCOMMAND.PIPELINE:  # load the configuration file
        config = MavisConfig.read(args.config)
        config.output = args.output
        config.skip_stage = args.skip_stage
        rargs = config.reference
        args = config

    # set all reference files to their absolute paths to make tracking them down later easier
    for arg in ['output', 'reference_genome', 'template_metadata', 'annotations', 'masking', 'aligner_reference',
                'dgv_annotation']:
        try:
            rargs[arg] = os.path.abspath(rargs[arg])
            if arg != 'output' and not os.path.isfile(rargs[arg]):
                raise OSError('input reference file does not exist', arg, rargs[arg])
        except AttributeError:
            pass

    # try checking the input files exist
    try:
        for fname in args.inputs:
            if len(bash_expands(fname)) < 1:
                raise OSError('input file does not exist', fname)
    except AttributeError:
        pass
    try:
        if len(bash_expands(args.input)) < 1:
            raise OSError('input file does not exist', args.input)
    except AttributeError:
        pass

    # load the reference files if they have been given and reset the arguments to hold the original file name and the
    # loaded data
    try:
        if any([
            pstep == SUBCOMMAND.CLUSTER and args.uninformative_filter,
            pstep == SUBCOMMAND.PIPELINE and config.cluster.uninformative_filter,
            pstep == SUBCOMMAND.VALIDATE and args.protocol == PROTOCOL.TRANS,
            pstep == SUBCOMMAND.PIPELINE and config.has_transcriptome() and SUBCOMMAND.VALIDATE not in config.skip_stage,
            pstep == SUBCOMMAND.PAIR or pstep == SUBCOMMAND.ANNOTATE or pstep == SUBCOMMAND.SUMMARY,
            pstep == SUBCOMMAND.OVERLAY
        ]):
            log('loading:', rargs.annotations)
            rargs.annotations_filename = rargs.annotations
            rargs.annotations = load_annotations(rargs.annotations)
        else:
            rargs.annotations_filename = rargs.annotations
            rargs.annotations = None
    except AttributeError as err:
        pass
    # reference genome
    try:
        if pstep in [SUBCOMMAND.VALIDATE, SUBCOMMAND.ANNOTATE]:
            log('loading:', rargs.reference_genome)
            rargs.reference_genome_filename = rargs.reference_genome
            rargs.reference_genome = load_reference_genome(rargs.reference_genome)
        else:
            rargs.reference_genome_filename = rargs.reference_genome
            rargs.reference_genome = None
    except AttributeError:
        pass

    # masking file
    try:
        if pstep in [SUBCOMMAND.VALIDATE, SUBCOMMAND.CLUSTER, SUBCOMMAND.PIPELINE]:
            log('loading:', rargs.masking)
            rargs.masking_filename = rargs.masking
            rargs.masking = load_masking_regions(rargs.masking)
        else:
            rargs.masking_filename = rargs.masking
            rargs.masking = None
    except AttributeError:
        pass

    # dgv annotation
    try:
        if pstep == SUBCOMMAND.SUMMARY:
            log('loading:', rargs.dgv_annotation)
            rargs.dgv_annotation_filename = rargs.dgv_annotation
            rargs.dgv_annotation = load_masking_regions(rargs.dgv_annotation)
        else:
            rargs.dgv_annotation_filename = rargs.dgv_annotation
            rargs.dgv_annotation = None
    except AttributeError:
        pass

    # template metadata
    try:
        if pstep == SUBCOMMAND.ANNOTATE:
            log('loading:', rargs.template_metadata)
            rargs.template_metadata_filename = rargs.template_metadata
            rargs.template_metadata = load_templates(rargs.template_metadata)
        else:
            rargs.template_metadata_filename = rargs.template_metadata
            rargs.template_metadata = None
    except AttributeError:
        pass
    # decide which main function to execute
    if pstep == SUBCOMMAND.CLUSTER:
        cluster_main(**args, start_time=start_time)
    elif pstep == SUBCOMMAND.VALIDATE:
        validate_main(**args, start_time=start_time)
    elif pstep == SUBCOMMAND.ANNOTATE:
        annotate_main(**args, start_time=start_time)
    elif pstep == SUBCOMMAND.PAIR:
        pairing_main(**args, start_time=start_time)
    elif pstep == SUBCOMMAND.SUMMARY:
        summary_main(**args, start_time=start_time)
    elif pstep == SUBCOMMAND.CONVERT:
        convert_main(**args)
    elif pstep == SUBCOMMAND.OVERLAY:
        overlay_main(**args)
    else:  # PIPELINE
        main_pipeline(args)

    duration = int(time.time()) - start_time
    hours = duration - duration % 3600
    minutes = duration - hours - (duration - hours) % 60
    seconds = duration - hours - minutes
    log(
        'run time (hh/mm/ss): {}:{:02d}:{:02d}'.format(hours // 3600, minutes // 60, seconds),
        time_stamp=False)
    log('run time (s): {}'.format(duration), time_stamp=False)
    return EXIT_OK


if __name__ == '__main__':
    exit(main())
