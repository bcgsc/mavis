#!python
import argparse
import os
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
from .config import augment_parser, MavisConfig, generate_config, CustomHelpFormatter, CONVERT_OPTIONS
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
    # read the config
    # set up the directory structure and run mavis
    batch_id = 'batch-' + str(uuid.uuid4())
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
                scriptname = script.write(os.path.join(outputdir, 'submit.sh'))

                submitall.append('vjob{}=$({} {})'.format(jobid_var_index, SCHEDULER[config.schedule.scheduler].submit, scriptname))
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
            script = SubmissionScript(command, config.schedule.scheduler, **options)
            scriptname = script.write(os.path.join(outputdir, 'submit.sh'))
            prevjob = '${{vjob{}##* }}'.format(jobid_var_index)
            submitall.append('ajob{}=$({} {} {})'.format(
                jobid_var_index, SCHEDULER[config.schedule.scheduler].submit,
                SCHEDULER[config.schedule.scheduler].dependency(prevjob),
                scriptname))
            outputfile = os.path.join(outputdir, ANNOTATION_PASS_PATTERN)
            pairing_inputs.append(outputfile)
            job_name_by_output[outputfile] = options['jobname']
            jobid_var_index += 1

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
    script = SubmissionScript(command, config.schedule.scheduler, **options)
    scriptname = script.write(os.path.join(outputdir, 'submit.sh'))

    submitall.append('jobid=$({} {} {})'.format(
        SCHEDULER[config.schedule.scheduler].submit,
        SCHEDULER[config.schedule.scheduler].dependency(
            ':'.join(['${{ajob{}##* }}'.format(i) for i in range(0, jobid_var_index)])),
        scriptname))

    # set up scripts for the summary held on the pairing job
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
    script = SubmissionScript(command, config.schedule.scheduler, **options)
    scriptname = script.write(os.path.join(outputdir, 'submit.sh'))

    submitall.append('{} {} {}'.format(
        SCHEDULER[config.schedule.scheduler].submit,
        SCHEDULER[config.schedule.scheduler].dependency('${jobid##* }'),
        scriptname))

    # now write a script at the top level to submit all
    submitallfile = os.path.join(config.output, 'submit_pipeline_{}.sh'.format(batch_id))
    log('writing:', submitallfile)
    with open(submitallfile, 'w') as fh:
        for line in submitall:
            fh.write(line + '\n')


def convert_main(inputs, output, file_type, strand_specific=False, assume_no_untemplated=True):
    bpp_results = []
    for filename in inputs:
        bpp_results.extend(convert_tool_output(filename, file_type, strand_specific, log, True, assume_no_untemplated=assume_no_untemplated))
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
            augment_parser(['strand_specific', 'assume_no_untemplated'], optional)

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
        convert_main(args.inputs, args.output, args.file_type, args.strand_specific, assume_no_untemplated=args.assume_no_untemplated)
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
