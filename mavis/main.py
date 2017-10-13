#!python
import argparse
from datetime import datetime
import math
import os
import random
import re
import subprocess
import sys
import time

import tab
import uuid

from . import __version__
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .annotate.file_io import load_annotations, load_masking_regions, load_reference_genome, load_templates
from .annotate.main import main as annotate_main
from .bam.read import get_samtools_version
from .blat import get_blat_version
from .checker import check_completion
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .cluster.main import main as cluster_main
from .cluster.main import split_clusters
from .config import augment_parser, get_metavar, LibraryConfig, MavisConfig, write_config
from .constants import PIPELINE_STEP, PROTOCOL
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .pairing.main import main as pairing_main
from .submit import SubmissionScript, SCHEDULER
from .submit import OPTIONS as SUBMIT_OPTIONS
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .summary.main import main as summary_main
from .tools import convert_tool_output, SUPPORTED_TOOL
from .util import bash_expands, log, log_arguments, MavisNamespace, mkdirp, output_tabbed_file, get_env_variable
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
    args.update(config.validation.items())
    args.update({k: v for k, v in libconf.items() if k in args})

    command = ['{} {}'.format(PROGNAME, PIPELINE_STEP.VALIDATE)]
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
        'min_orf_size': config.annotation.min_orf_size,
        'max_orf_cap': config.annotation.max_orf_cap,
        'library': libconf.library,
        'protocol': libconf.protocol,
        'min_domain_mapping_match': config.annotation.min_domain_mapping_match,
        'domain_name_regex_filter': config.illustrate.domain_name_regex_filter,
        'max_proximity': config.cluster.max_proximity
    }
    args.update({k: v for k, v in libconf.items() if k in args})
    command = ['{} {}'.format(PROGNAME, PIPELINE_STEP.ANNOTATE)]
    for argname, value in args.items():
        if isinstance(value, str):
            command.append('--{} "{}"'.format(argname, value))
        else:
            command.append('--{} {}'.format(argname, value))
    command.append('--inputs {}'.format(inputfile))
    command.append('--output {}'.format(outputdir))
    return ' \\\n\t'.join(command) + '\n'


def run_conversion(config, libconf, conversion_dir):
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
                    output_tabbed_file(convert_tool_output(*command[1:], log=log), output_filename)
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
    # read the config
    # set up the directory structure and run mavis
    batch_id = 'batch-' + str(uuid.uuid4())
    conversion_dir = mkdirp(os.path.join(config.output, 'converted_inputs'))
    pairing_inputs = []
    submission_scripts = []

    def get_prefix(filename):
        return re.sub(r'\.[^\.]+$', '', os.path.basename(filename))

    job_name_by_output = {}

    for libconf in config.libraries.values():
        base = os.path.join(config.output, '{}_{}_{}'.format(libconf.library, libconf.disease_status, libconf.protocol))
        log('setting up the directory structure for', libconf.library, 'as', base)
        libconf.inputs = run_conversion(config, libconf, conversion_dir)

        # run the cluster stage
        cluster_output = mkdirp(os.path.join(base, PIPELINE_STEP.CLUSTER))  # creates the output dir
        merge_args = {'batch_id': batch_id, 'output': cluster_output}
        merge_args['split_only'] = PIPELINE_STEP.CLUSTER in config.skip_stage
        merge_args.update(config.reference.items())
        merge_args.update(config.cluster.items())
        merge_args.update(libconf.items())
        log('clustering', '(split only)' if merge_args['split_only'] else '')
        inputs = cluster_main(log_args=True, **merge_args)

        for inputfile in inputs:
            prefix = get_prefix(inputfile)  # will be batch id + job number

            if PIPELINE_STEP.VALIDATE not in config.skip_stage:
                outputdir = mkdirp(os.path.join(base, PIPELINE_STEP.VALIDATE, prefix))
                command = build_validate_command(config, libconf, inputfile, outputdir)
                # build the submission script
                options = {k: config.schedule[k] for k in SUBMIT_OPTIONS}
                options['stdout'] = outputdir
                options['jobname'] = 'MV_{}_{}'.format(libconf.library, prefix)

                if libconf.is_trans():
                    options['memory_limit'] = config.schedule.trans_validation_memory
                else:
                    options['memory_limit'] = config.schedule.validation_memory
                script = SubmissionScript(command, config.schedule.scheduler, **options)
                submission_scripts.append(script.write(os.path.join(outputdir, 'submit.sh')))

                # for setting up subsequent jobs and holds
                outputfile = os.path.join(outputdir, VALIDATION_PASS_PATTERN)
                job_name_by_output[outputfile] = options['jobname']
                inputfile = outputfile

            # annotation cannot be skipped
            outputdir = mkdirp(os.path.join(base, PIPELINE_STEP.ANNOTATE, prefix))
            command = build_annotate_command(config, libconf, inputfile, outputdir)

            options = {k: config.schedule[k] for k in SUBMIT_OPTIONS}
            options['stdout'] = outputdir
            options['jobname'] = 'MA_{}_{}'.format(libconf.library, prefix)
            options['memory_limit'] = config.schedule.annotation_memory
            if inputfile in job_name_by_output:
                options['dependency'] = job_name_by_output[inputfile]
            script = SubmissionScript(command, config.schedule.scheduler, **options)
            submission_scripts.append(script.write(os.path.join(outputdir, 'submit.sh')))

            outputfile = os.path.join(outputdir, ANNOTATION_PASS_PATTERN)
            pairing_inputs.append(outputfile)
            job_name_by_output[outputfile] = options['jobname']

    # set up scripts for the pairing held on all of the annotation jobs
    outputdir = mkdirp(os.path.join(config.output, PIPELINE_STEP.PAIR))
    args = config.pairing.flatten()
    args.update({
        'output': outputdir,
        'annotations': config.reference.annotations_filename
    })
    command = ['{} {}'.format(PROGNAME, PIPELINE_STEP.PAIR)]
    for arg, value in sorted(args.items()):
        if isinstance(value, str):
            command.append('--{} "{}"'.format(arg, value))
        else:
            command.append('--{} {}'.format(arg, value))
    command.append('--inputs {}'.format(' \\\n\t'.join(pairing_inputs)))
    command = ' \\\n\t'.join(command)

    options = {k: config.schedule[k] for k in SUBMIT_OPTIONS}
    options['stdout'] = outputdir
    options['jobname'] = 'MP_{}'.format(batch_id)
    options['dependency'] = ','.join(job_name_by_output[o] for o in pairing_inputs)
    script = SubmissionScript(command, config.schedule.scheduler, **options)
    submission_scripts.append(script.write(os.path.join(outputdir, 'submit.sh')))

    # set up scripts for the summary held on the pairing job
    pairing_jobname = script.jobname
    outputdir = mkdirp(os.path.join(config.output, PIPELINE_STEP.SUMMARY))
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
    command = ['{} {}'.format(PROGNAME, PIPELINE_STEP.SUMMARY)]
    for arg, value in sorted(args.items()):
        if isinstance(value, str):
            command.append('--{} "{}"'.format(arg, value))
        else:
            command.append('--{} {}'.format(arg, value))
    command = ' \\\n\t'.join(command)

    options = {k: config.schedule[k] for k in SUBMIT_OPTIONS}
    options['stdout'] = outputdir
    options['jobname'] = 'MS_{}'.format(batch_id)
    options['dependency'] = pairing_jobname
    script = SubmissionScript(command, config.schedule.scheduler, **options)
    submission_scripts.append(script.write(os.path.join(outputdir, 'submit.sh')))

    # now write a script at the top level to submit all
    submitall = os.path.join(config.output, 'submit_pipeline_{}.sh'.format(batch_id))
    log('writing:', submitall)
    with open(submitall, 'w') as fh:
        for script in submission_scripts:
            fh.write('{} {}\n'.format(SCHEDULER[config.schedule.scheduler].submit, script))


def generate_config(parser, required, optional):
    """
    Args:
        parser (argparse.ArgumentParser): the main parser
        required: the argparse required arguments group
        optional: the argparse optional arguments group
    """
    # the config sub  program is used for writing pipeline configuration files
    required.add_argument('-w', '--write', help='path to the new configuration file', required=True, metavar='FILEPATH')
    optional.add_argument(
        '--library', nargs=5,
        metavar=('<name>', '(genome|transcriptome)', '<diseased|normal>', '</path/to/bam/file>', '<strand_specific>'),
        action='append', help='configuration for libraries to be analyzed by mavis', default=[])
    optional.add_argument(
        '--input', help='path to an input file or filter for mavis followed by the library names it '
        'should be used for', nargs='+', action='append', default=[], metavar='FILEPATH'
    )
    optional.add_argument(
        '--assign', help='library name followed by path(s) to input file(s) or filter names. This represents the list'
        ' of inputs that should be used for the library', nargs='+', default=[], action='append')
    optional.add_argument(
        '--best_transcripts_only', default=get_env_variable('best_transcripts_only', True), metavar=get_metavar(bool),
        type=tab.cast_boolean, help='compute from best transcript models only')
    optional.add_argument(
        '--genome_bins', default=get_env_variable('genome_bins', 100), type=int, metavar=get_metavar(int),
        help='number of bins/samples to use in calculating the fragment size stats for genomes')
    optional.add_argument(
        '--transcriptome_bins', default=get_env_variable('transcriptome_bins', 5000), type=int, metavar=get_metavar(int),
        help='number of genes to use in calculating the fragment size stats for genomes')
    optional.add_argument(
        '--distribution_fraction', default=get_env_variable('distribution_fraction', 0.97), type=float, metavar=get_metavar(float),
        help='the proportion of the distribution of calculated fragment sizes to use in determining the stdev')
    optional.add_argument(
        '--verbose', default=get_env_variable('verbose', False), type=tab.cast_boolean, metavar=get_metavar(bool),
        help='verbosely output logging information')
    optional.add_argument(
        '--convert', nargs=4, default=[],
        metavar=('<alias>', '</path/to/input/file>', '({})'.format('|'.join(SUPPORTED_TOOL.values())), '<stranded>'),
        help='input file conversion for internally supported tools', action='append')
    optional.add_argument(
        '--external_conversion', metavar=('<alias>', '<"command">'), nargs=2, default=[],
        help='alias for use in inputs and full command (quoted)', action='append')
    optional.add_argument(
        '--no_defaults', default=False, action='store_true', help='do not write current defaults to the config output')
    augment_parser(['annotations'], required, optional)
    args = parser.parse_args()
    if args.distribution_fraction < 0 or args.distribution_fraction > 1:
        raise ValueError('distribution_fraction must be a value between 0-1')
    log('MAVIS: {}'.format(__version__))
    log_arguments(args.__dict__)

    # process the libraries by input argument (--input)
    inputs_by_lib = {k[0]: set() for k in args.library}
    for arg_list in args.input:
        if len(arg_list) < 2:
            raise ValueError('--input requires 2+ arguments', arg_list)
        inputfile = arg_list[0]
        for lib in arg_list[1:]:
            if lib not in inputs_by_lib:
                raise KeyError(
                    '--input specified a library that was not configured. Please input all libraries using '
                    'the --library flag', lib)
            inputs_by_lib[lib].add(inputfile)
    # process the inputs by library argument (--assign)
    for arg_list in args.assign:
        if len(arg_list) < 2:
            raise ValueError('--assign requires 2+ arguments', arg_list)
        lib = arg_list[0]
        if lib not in inputs_by_lib:
            raise KeyError(
                '--assign specified a library that was not configured. Please input all libraries using '
                'the --library flag', lib)
        inputs_by_lib[lib].update(arg_list[1:])

    libs = []
    # load the annotations if we need them
    if any([p == 'transcriptome' for l, p, d, b, s in args.library]):
        log('loading the reference annotations file', args.annotations)
        args.annotations_filename = args.annotations
        args.annotations = load_annotations(args.annotations, best_transcripts_only=args.best_transcripts_only)

    for lib, protocol, diseased, bam, stranded in args.library:
        if lib not in inputs_by_lib:
            raise KeyError('not input was given for the library', lib)
        log('generating the config section for:', lib)
        temp = LibraryConfig.build(
            library=lib, protocol=protocol, bam_file=bam, inputs=inputs_by_lib[lib], strand_specific=stranded,
            disease_status=diseased, annotations=args.annotations, log=log,
            sample_size=args.genome_bins if protocol == PROTOCOL.GENOME else args.transcriptome_bins,
            distribution_fraction=args.distribution_fraction
        )
        libs.append(temp)
    convert = {}
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

    for alias, inputfile, toolname, stranded in args.convert:
        if alias in convert:
            raise KeyError('duplicate alias names are not allowed', alias)
        stranded = str(tab.cast_boolean(stranded))
        SUPPORTED_TOOL.enforce(toolname)
        convert[alias] = ['convert_tool_output', inputfile, toolname, stranded]
    write_config(args.write, include_defaults=not args.no_defaults, libraries=libs, conversions=convert, log=log)


def convert_main(inputs, output, file_type, strand_specific=False):
    bpp_results = []
    for filename in inputs:
        bpp_results.extend(convert_tool_output(filename, file_type, strand_specific, log, True))
    output_filename = 'mavis_{}_converted.tab'.format(file_type)
    output_tabbed_file(bpp_results, os.path.join(output, output_filename))


def main():
    def usage(err=None, detail=False):
        umsg = '\nusage: {} {{{}}} [-h] [-v]'.format(PROGNAME, ','.join(sorted(PIPELINE_STEP.values())))
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

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    augment_parser(['help', 'version'], optional)

    if pstep == PIPELINE_STEP.CONFIG:
        generate_config(parser, required, optional)
        return EXIT_OK
    else:
        required.add_argument('-o', '--output', help='path to the output directory', required=True)

        if pstep == PIPELINE_STEP.PIPELINE:
            required.add_argument('config', help='path to the input pipeline configuration file', metavar='FILEPATH')
            optional.add_argument(
                '--skip_stage', choices=[PIPELINE_STEP.CLUSTER, PIPELINE_STEP.VALIDATE], action='append', default=[],
                help='Use flag once per stage to skip. Can skip clustering or validation or both')

        elif pstep == PIPELINE_STEP.CONVERT:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            required.add_argument(
                '--file_type', choices=sorted([t for t in SUPPORTED_TOOL.values() if t != 'mavis']),
                required=True, help='Indicates the input file type to be parsed')
            augment_parser(['strand_specific'], optional)

        elif pstep == PIPELINE_STEP.CLUSTER:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(['library', 'protocol', 'strand_specific', 'disease_status', 'annotations', 'masking'], required, optional)
            augment_parser(CLUSTER_DEFAULTS.keys(), optional)

        elif pstep == PIPELINE_STEP.VALIDATE:
            required.add_argument('-n', '--input', help='path to the input file', required=True, metavar='FILEPATH')
            augment_parser(
                ['library', 'protocol', 'bam_file', 'read_length', 'stdev_fragment_size', 'median_fragment_size'] +
                ['strand_specific', 'annotations', 'reference_genome', 'aligner_reference', 'masking'],
                required, optional
            )
            augment_parser(VALIDATION_DEFAULTS.keys(), optional)

        elif pstep == PIPELINE_STEP.ANNOTATE:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(
                ['library', 'protocol', 'annotations', 'reference_genome', 'masking', 'template_metadata'],
                required, optional
            )
            augment_parser(['max_proximity'], optional)
            augment_parser(list(ANNOTATION_DEFAULTS.keys()) + list(ILLUSTRATION_DEFAULTS.keys()), optional)

        elif pstep == PIPELINE_STEP.PAIR:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            optional.add_argument(
                '-f', '--product_sequence_files', nargs='+', help='paths to fasta files with product sequences', metavar='FILEPATH',
                required=False, default=[])
            augment_parser(['annotations'], required, optional)
            augment_parser(['max_proximity'] + list(PAIRING_DEFAULTS.keys()), optional)

        elif pstep == PIPELINE_STEP.SUMMARY:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True, metavar='FILEPATH')
            augment_parser(
                ['annotations', 'dgv_annotation', 'flanking_call_distance', 'split_call_distance', 'contig_call_distance', 'spanning_call_distance'],
                required, optional
            )
            augment_parser(SUMMARY_DEFAULTS.keys(), optional)
        elif pstep == PIPELINE_STEP.CHECKER:

            args = parser.parse_args()
            success_flag = check_completion(args.output)
            return EXIT_OK if success_flag else EXIT_ERROR

        else:
            raise NotImplementedError('invalid value for <pipeline step>', pstep)
    args = MavisNamespace(**parser.parse_args().__dict__)
    args.samtools_version = get_samtools_version()
    args.blat_version = get_blat_version()

    log('MAVIS: {}'.format(__version__))
    log_arguments(args.__dict__)
    rargs = args

    if pstep == PIPELINE_STEP.PIPELINE:  # load the configuration file
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
            pstep == PIPELINE_STEP.CLUSTER and args.uninformative_filter,
            pstep == PIPELINE_STEP.PIPELINE and config.cluster.uninformative_filter,
            pstep == PIPELINE_STEP.VALIDATE and args.protocol == PROTOCOL.TRANS,
            pstep == PIPELINE_STEP.PIPELINE and config.has_transcriptome(),
            pstep == PIPELINE_STEP.PAIR or pstep == PIPELINE_STEP.ANNOTATE or pstep == PIPELINE_STEP.SUMMARY
        ]):
            log('loading:', rargs.annotations)
            rargs.annotations_filename = rargs.annotations
            rargs.annotations = load_annotations(rargs.annotations)
        else:
            rargs.annotations_filename = rargs.annotations
            rargs.annotations = None
    except AttributeError:
        pass
    # reference genome
    try:
        if pstep in [PIPELINE_STEP.VALIDATE, PIPELINE_STEP.ANNOTATE]:
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
        if pstep in [PIPELINE_STEP.VALIDATE, PIPELINE_STEP.CLUSTER, PIPELINE_STEP.PIPELINE]:
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
        if pstep == PIPELINE_STEP.SUMMARY:
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
        if pstep == PIPELINE_STEP.ANNOTATE:
            log('loading:', rargs.template_metadata)
            rargs.template_metadata_filename = rargs.template_metadata
            rargs.template_metadata = load_templates(rargs.template_metadata)
        else:
            rargs.template_metadata_filename = rargs.template_metadata
            rargs.template_metadata = None
    except AttributeError:
        pass
    # decide which main function to execute
    if pstep == PIPELINE_STEP.CLUSTER:
        cluster_main(**args)
    elif pstep == PIPELINE_STEP.VALIDATE:
        validate_main(**args)
    elif pstep == PIPELINE_STEP.ANNOTATE:
        annotate_main(**args)
    elif pstep == PIPELINE_STEP.PAIR:
        pairing_main(**args)
    elif pstep == PIPELINE_STEP.SUMMARY:
        summary_main(**args)
    elif pstep == PIPELINE_STEP.CONVERT:
        convert_main(args.inputs, args.output, args.file_type, args.strand_specific)
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
