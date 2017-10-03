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

from . import __version__
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .annotate.file_io import load_annotations, load_masking_regions, load_reference_genome, load_templates
from .annotate.main import main as annotate_main
from .bam.read import get_samtools_version
from .blat import get_blat_version
from .checker import check_completion
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .cluster.main import main as cluster_main
from .config import augment_parser, get_env_variable, LibraryConfig, MavisConfig, write_config
from .constants import PIPELINE_STEP, PROTOCOL
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .pairing.main import main as pairing_main
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from .summary.main import main as summary_main
from .tools import convert_tool_output, SUPPORTED_TOOL
from .util import bash_expands, log, log_arguments, MavisNamespace, mkdirp, output_tabbed_file
from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from .validate.main import main as validate_main


VALIDATION_PASS_SUFFIX = '.validation-passed.tab'
PROGNAME = 'mavis'
JOB_TASK = "$SGE_TASK_ID"
EXIT_OK = 0
EXIT_ERROR = 1

QSUB_HEADER = """#!/bin/bash
#$ -V
#$ -N {name}
#$ -q {queue}
#$ -o {output}
#$ -l mem_free={memory}G,mem_token={memory}G,h_vmem={memory}G
#$ -j y"""


def main_pipeline(config):
    # read the config
    # set up the directory structure and run mavis
    annotation_jobs = []
    rand = int(random.random() * math.pow(10, 10))
    conversion_dir = mkdirp(os.path.join(config.output, 'converted_inputs'))
    pairing_inputs = []
    for libconf in config.libraries.values():
        base = os.path.join(config.output, '{}_{}_{}'.format(libconf.library, libconf.disease_status, libconf.protocol))
        log('setting up the directory structure for', libconf.library, 'as', base)
        cluster_output = mkdirp(os.path.join(base, PIPELINE_STEP.CLUSTER))
        validation_output = mkdirp(os.path.join(base, PIPELINE_STEP.VALIDATE))
        annotation_output = mkdirp(os.path.join(base, PIPELINE_STEP.ANNOTATE))
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
        libconf.inputs = inputs
        # run the merge
        log('clustering')
        merge_args = {}
        merge_args.update(config.reference.items())
        merge_args.update(config.cluster.items())
        merge_args.update(libconf.items())
        merge_args['output'] = cluster_output
        output_files = cluster_main(log_args=True, **merge_args)
        job_task = '1' if len(output_files) == 1 else JOB_TASK
        if not output_files:
            log('warning: no inputs after clustering. Will not set up other pipeline steps')
            continue
        merge_file_prefix = None
        for filename in output_files:
            m = re.match(r'^(?P<prefix>.*\D)\d+.tab$', filename)
            if not m:
                raise UserWarning('cluster file did not match expected format', filename)
            if merge_file_prefix is None:
                merge_file_prefix = m.group('prefix')
            elif merge_file_prefix != m.group('prefix'):
                raise UserWarning('merge file prefixes are not consistent', output_files)

        # now set up the qsub script for the validation and the held job for the annotation
        validation_args = {
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
            'stranded_bam': libconf.stranded_bam
        }
        validation_args.update(config.validation.items())
        validation_args.update({k: v for k, v in libconf.items() if k in validation_args})

        qsub = os.path.join(validation_output, 'qsub.sh')
        validation_jobname = '{}_{}_{}_{}'.format(PIPELINE_STEP.VALIDATE, libconf.library, libconf.protocol, rand)

        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=config.schedule.queue,
                    memory=config.schedule.trans_validation_memory if libconf.is_trans() else config.schedule.validation_memory,
                    name=validation_jobname,
                    output=validation_output
                ) + '\n')
            if job_task != '1':
                fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            temp = [
                '--{} {}'.format(k, v) for k, v in validation_args.items() if not isinstance(v, str) and v is not None]
            temp.extend(
                ['--{} "{}"'.format(k, v) for k, v in validation_args.items() if isinstance(v, str) and v is not None])
            validation_args = temp
            validation_args.append('-n {}{}.tab'.format(merge_file_prefix, job_task))
            fh.write('{} {} {}'.format(PROGNAME, PIPELINE_STEP.VALIDATE, ' \\\n\t'.join(validation_args)))
            fh.write(
                ' \\\n\t--output {}\n'.format(
                    os.path.join(validation_output, os.path.basename(merge_file_prefix) + job_task)))

        # set up the annotations job
        # for all files with the right suffix
        annotation_args = {
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
        annotation_args.update({k: v for k, v in libconf.items() if k in annotation_args})
        temp = [
            '--{} {}'.format(k, v) for k, v in annotation_args.items() if not isinstance(v, str) and v is not None]
        temp.extend(
            ['--{} "{}"'.format(k, v) for k, v in annotation_args.items() if isinstance(v, str) and v is not None])
        annotation_args = temp
        annotation_args.append('--inputs {}/{}{}/*{}'.format(
            validation_output, os.path.basename(merge_file_prefix), job_task, VALIDATION_PASS_SUFFIX))
        qsub = os.path.join(annotation_output, 'qsub.sh')
        annotation_jobname = '{}_{}_{}_{}'.format(PIPELINE_STEP.ANNOTATE, libconf.library, libconf.protocol, rand)
        annotation_jobs.append(annotation_jobname)
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=config.schedule.queue,
                    memory=libconf.get('annotation_memory', config.schedule.annotation_memory),
                    name=annotation_jobname,
                    output=annotation_output
                ) + '\n')
            if job_task != '1':
                fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            fh.write('#$ -hold_jid {}\n'.format(validation_jobname))
            fh.write('{} {} {}'.format(PROGNAME, PIPELINE_STEP.ANNOTATE, ' \\\n\t'.join(annotation_args)))
            fh.write(
                ' \\\n\t--output {}\n'.format(
                    os.path.join(annotation_output, os.path.basename(merge_file_prefix) + job_task)))
            for sge_task_id in range(1, len(output_files) + 1):
                pairing_inputs.append(os.path.join(
                    annotation_output, os.path.basename(merge_file_prefix) + str(sge_task_id), 'annotations.tab'))

    # set up scripts for the pairing held on all of the annotation jobs
    pairing_output = mkdirp(os.path.join(config.output, PIPELINE_STEP.PAIR))
    pairing_args = config.pairing.flatten()
    pairing_args.update({
        'output': pairing_output,
        'annotations': config.reference.annotations_filename
    })
    temp = ['--{} {}'.format(k, v) for k, v in pairing_args.items() if not isinstance(v, str) and v is not None]
    temp.extend(['--{} "{}"'.format(k, v) for k, v in pairing_args.items() if isinstance(v, str) and v is not None])
    temp.append('--inputs {}'.format(' \\\n\t'.join(pairing_inputs)))
    pairing_args = temp
    qsub = os.path.join(pairing_output, 'qsub.sh')
    pairing_jobname = 'mavis_{}_{}'.format(PIPELINE_STEP.PAIR, rand)
    with open(qsub, 'w') as fh:
        log('writing:', qsub)
        fh.write(
            QSUB_HEADER.format(
                queue=config.schedule.queue, memory=config.schedule.memory, name=pairing_jobname, output=pairing_output
            ) + '\n')
        fh.write('#$ -hold_jid {}\n'.format(','.join(annotation_jobs)))
        fh.write('{} {} {}\n'.format(PROGNAME, PIPELINE_STEP.PAIR, ' \\\n\t'.join(pairing_args)))

    # set up scripts for the summary held on the pairing job
    summary_output = mkdirp(os.path.join(config.output, PIPELINE_STEP.SUMMARY))
    summary_args = dict(
        output=summary_output,
        flanking_call_distance=config.pairing.flanking_call_distance,
        split_call_distance=config.pairing.split_call_distance,
        contig_call_distance=config.pairing.contig_call_distance,
        spanning_call_distance=config.pairing.spanning_call_distance,
        dgv_annotation=config.reference.dgv_annotation_filename,
        annotations=config.reference.annotations_filename
    )
    summary_args.update(config.summary.items())
    temp = ['--{} {}'.format(k, v) for k, v in summary_args.items() if not isinstance(v, str) and v is not None]
    temp.extend(['--{} "{}"'.format(k, v) for k, v in summary_args.items() if isinstance(v, str) and v is not None])
    temp.append('--inputs {}'.format(os.path.join(config.output, 'pairing/mavis_paired*.tab')))
    summary_args = temp
    qsub = os.path.join(summary_output, 'qsub.sh')
    with open(qsub, 'w') as fh:
        log('writing:', qsub)
        fh.write(
            QSUB_HEADER.format(
                queue=config.schedule.queue,
                memory=config.schedule.memory,
                name='mavis_{}'.format(PIPELINE_STEP.SUMMARY),
                output=summary_output
            ) + '\n')
        fh.write('#$ -hold_jid {}\n'.format(pairing_jobname))
        fh.write('{} {} {}\n'.format(PROGNAME, PIPELINE_STEP.SUMMARY, ' \\\n\t'.join(summary_args)))


def generate_config(parser, required, optional):
    """
    Args:
        parser (argparse.ArgumentParser): the main parser
        required: the argparse required arguments group
        optional: the argparse optional arguments group
    """
    # the config sub  program is used for writing pipeline configuration files
    required.add_argument('-w', '--write', help='path to the new configuration file', required=True)
    optional.add_argument(
        '--library', nargs=5,
        metavar=('<name>', '(genome|transcriptome)', '<diseased|normal>', '</path/to/bam/file>', '<stranded_bam>'),
        action='append', help='configuration for libraries to be analyzed by mavis', default=[])
    optional.add_argument(
        '--input', help='path to an input file or filter for mavis followed by the library names it '
        'should be used for', nargs='+', action='append', default=[]
    )
    optional.add_argument(
        '--assign', help='library name followed by path(s) to input file(s) or filter names. This represents the list'
        ' of inputs that should be used for the library', nargs='+', default=[], action='append')
    optional.add_argument(
        '--best_transcripts_only', default=get_env_variable('best_transcripts_only', True),
        type=tab.cast_boolean, help='compute from best transcript models only')
    optional.add_argument(
        '--genome_bins', default=get_env_variable('genome_bins', 100), type=int,
        help='number of bins/samples to use in calculating the fragment size stats for genomes')
    optional.add_argument(
        '--transcriptome_bins', default=get_env_variable('transcriptome_bins', 5000), type=int,
        help='number of genes to use in calculating the fragment size stats for genomes')
    optional.add_argument(
        '--distribution_fraction', default=get_env_variable('distribution_fraction', 0.97), type=float,
        help='the proportion of the distribution of calculated fragment sizes to use in determining the stdev')
    optional.add_argument(
        '--verbose', default=get_env_variable('verbose', False), type=tab.cast_boolean,
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
            library=lib, protocol=protocol, bam_file=bam, inputs=inputs_by_lib[lib], stranded_bam=stranded,
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


def time_diff(start, end):
    """
    >>> time_diff('2017-04-02 15:10:48.607195', '2017-04-03 17:00:32.671809')
    25.83
    >>> time_diff('2017-04-03 15:10:48.607195', '2017-04-03 17:00:32.671809')
    1.83
    """
    start_time = datetime.strptime(start, '%Y-%m-%d %H:%M:%S.%f')
    end_time = datetime.strptime(end, '%Y-%m-%d %H:%M:%S.%f')
    elapsed_time = end_time - start_time
    return elapsed_time.seconds / 3600


def main():
    def usage(err=None, detail=False):
        umsg = '\nusage: {} {{cluster,validate,annotate,pairing,summary,pipeline,config,checker}} [-h] [-v]'.format(PROGNAME)
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
            required.add_argument('config', help='path to the input pipeline configuration file')
        elif pstep == PIPELINE_STEP.CLUSTER:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            augment_parser(['library', 'protocol', 'stranded_bam', 'disease_status', 'annotations', 'masking'], required, optional)
            augment_parser(CLUSTER_DEFAULTS.keys(), optional)
        elif pstep == PIPELINE_STEP.VALIDATE:
            required.add_argument('-n', '--input', help='path to the input file', required=True)
            augment_parser(
                ['library', 'protocol', 'bam_file', 'read_length', 'stdev_fragment_size', 'median_fragment_size'] +
                ['stranded_bam', 'annotations', 'reference_genome', 'aligner_reference', 'masking'],
                required, optional
            )
            augment_parser(VALIDATION_DEFAULTS.keys(), optional)
        elif pstep == PIPELINE_STEP.ANNOTATE:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            augment_parser(
                ['library', 'protocol', 'annotations', 'reference_genome', 'masking', 'max_proximity', 'template_metadata'],
                required, optional
            )
            augment_parser(list(ANNOTATION_DEFAULTS.keys()) + list(ILLUSTRATION_DEFAULTS.keys()), optional)
        elif pstep == PIPELINE_STEP.PAIR:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            optional.add_argument(
                '-f', '--product_sequence_files', nargs='+', help='paths to fasta files with product sequences',
                required=False, default=[])
            augment_parser(['annotations'], required, optional)
            augment_parser(['max_proximity'] + list(PAIRING_DEFAULTS.keys()), optional)
        elif pstep == PIPELINE_STEP.SUMMARY:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
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

    # reference genome
    try:
        if pstep in [PIPELINE_STEP.VALIDATE, PIPELINE_STEP.ANNOTATE]:
            log('loading:', rargs.reference_genome)
            rargs.reference_genome_filename = rargs.reference_genome
            rargs.reference_genome = load_reference_genome(rargs.reference_genome)
        else:
            rargs.reference_genome_filename = rargs.reference_genome
            rargs.reference_genome = None
    except AttributeError as err:
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
