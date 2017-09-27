#!python
import os
import re
import sys
import random
import argparse
import TSV
import glob
import time
import subprocess
from .annotate import load_annotations, load_reference_genome, load_masking_regions, load_templates
from .validate.main import main as validate_main
from .cluster.main import main as cluster_main
from .pairing.main import main as pairing_main
from .annotate.main import main as annotate_main
from .summary.main import main as summary_main
from . import __version__

from .validate.constants import DEFAULTS as VALIDATION_DEFAULTS
from .cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from .annotate.constants import DEFAULTS as ANNOTATION_DEFAULTS
from .pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from .summary.constants import DEFAULTS as SUMMARY_DEFAULTS

from .config import augment_parser, write_config, LibraryConfig, MavisConfig, get_env_variable
from .util import log, mkdirp, log_arguments, output_tabbed_file, bash_expands, MavisNamespace
from .constants import PROTOCOL, PIPELINE_STEP, DISEASE_STATUS
from .bam.read import get_samtools_version
from .blat import get_blat_version
from .tools import convert_tool_output, SUPPORTED_TOOL
import math
from datetime import datetime


VALIDATION_PASS_SUFFIX = '.validation-passed.tab'
PROGNAME = 'mavis'

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
        cluster_output = mkdirp(os.path.join(base, 'clustering'))
        validation_output = mkdirp(os.path.join(base, 'validation'))
        annotation_output = mkdirp(os.path.join(base, 'annotation'))
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

        if len(output_files) == 0:
            log('warning: no inputs after clustering. Will not set up other pipeline steps')
            continue
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
        validation_jobname = 'validation_{}_{}_{}'.format(libconf.library, libconf.protocol, rand)

        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=config.schedule.queue,
                    memory=config.schedule.trans_validation_memory if libconf.is_trans() else config.schedule.validation_memory,
                    name=validation_jobname,
                    output=validation_output
                ) + '\n')
            fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            temp = [
                '--{} {}'.format(k, v) for k, v in validation_args.items() if not isinstance(v, str) and v is not None]
            temp.extend(
                ['--{} "{}"'.format(k, v) for k, v in validation_args.items() if isinstance(v, str) and v is not None])
            validation_args = temp
            validation_args.append('-n {}$SGE_TASK_ID.tab'.format(merge_file_prefix))
            fh.write('{} validate {}'.format(PROGNAME, ' \\\n\t'.join(validation_args)))
            fh.write(
                ' \\\n\t--output {}\n'.format(
                    os.path.join(validation_output, os.path.basename(merge_file_prefix) + '$SGE_TASK_ID')))

        # set up the annotations job
        # for all files with the right suffix
        annotation_args = {
            'reference_genome': config.reference.reference_genome_filename,
            'annotations': config.reference.annotations_filename,
            'template_metadata': config.reference.template_metadata_filename,
            'masking': config.reference.masking_filename,
            'min_orf_size': config.annotation.min_orf_size,
            'max_orf_cap': config.annotation.max_orf_cap,
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
        annotation_args.append('--inputs {}/{}$SGE_TASK_ID/*{}'.format(
            validation_output, os.path.basename(merge_file_prefix), VALIDATION_PASS_SUFFIX))
        qsub = os.path.join(annotation_output, 'qsub.sh')
        annotation_jobname = 'annotation_{}_{}_{}'.format(libconf.library, libconf.protocol, rand)
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
            fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            fh.write('#$ -hold_jid {}\n'.format(validation_jobname))
            fh.write('{} annotate {}'.format(PROGNAME, ' \\\n\t'.join(annotation_args)))
            fh.write(
                ' \\\n\t--output {}\n'.format(
                    os.path.join(annotation_output, os.path.basename(merge_file_prefix) + '$SGE_TASK_ID')))
            for sge_task_id in range(1, len(output_files) + 1):
                pairing_inputs.append(os.path.join(
                    annotation_output, os.path.basename(merge_file_prefix) + str(sge_task_id), 'annotations.tab'))

    # set up scripts for the pairing held on all of the annotation jobs
    pairing_output = mkdirp(os.path.join(config.output, 'pairing'))
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
    pairing_jobname = 'mavis_pairing_{}'.format(rand)
    with open(qsub, 'w') as fh:
        log('writing:', qsub)
        fh.write(
            QSUB_HEADER.format(
                queue=config.schedule.queue, memory=config.schedule.memory, name=pairing_jobname, output=pairing_output
            ) + '\n')
        fh.write('#$ -hold_jid {}\n'.format(','.join(annotation_jobs)))
        fh.write('{} pairing {}\n'.format(PROGNAME, ' \\\n\t'.join(pairing_args)))

    # set up scripts for the summary held on the pairing job
    summary_output = mkdirp(os.path.join(config.output, 'summary'))
    summary_args = dict(
        output=summary_output,
        filter_min_remapped_reads=config.summary.filter_min_remapped_reads,
        filter_min_spanning_reads=config.summary.filter_min_spanning_reads,
        filter_min_flanking_reads=config.summary.filter_min_flanking_reads,
        filter_min_flanking_only_reads=config.summary.filter_min_flanking_only_reads,
        filter_min_split_reads=config.summary.filter_min_split_reads,
        filter_min_linking_split_reads=config.summary.filter_min_linking_split_reads,
        flanking_call_distance=config.pairing.flanking_call_distance,
        split_call_distance=config.pairing.split_call_distance,
        contig_call_distance=config.pairing.contig_call_distance,
        spanning_call_distance=config.pairing.spanning_call_distance,
        dgv_annotation=config.reference.dgv_annotation_filename,
        annotations=config.reference.annotations_filename
    )
    temp = ['--{} {}'.format(k, v) for k, v in summary_args.items() if not isinstance(v, str) and v is not None]
    temp.extend(['--{} "{}"'.format(k, v) for k, v in summary_args.items() if isinstance(v, str) and v is not None])
    temp.append('--inputs {}'.format(os.path.join(config.output, 'pairing/mavis_paired*.tab')))
    summary_args = temp
    qsub = os.path.join(summary_output, 'qsub.sh')
    with open(qsub, 'w') as fh:
        log('writing:', qsub)
        fh.write(
            QSUB_HEADER.format(
                queue=config.schedule.queue, memory=config.schedule.memory, name='mavis_summary', output=summary_output
            ) + '\n')
        fh.write('#$ -hold_jid {}\n'.format(pairing_jobname))
        fh.write('{} summary {}\n'.format(PROGNAME, ' \\\n\t'.join(summary_args)))


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
        type=TSV.tsv_boolean, help='compute from best transcript models only')
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
        '--verbose', default=get_env_variable('verbose', False), type=TSV.tsv_boolean,
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
    augment_parser(required, optional, ['annotations'])
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
        l = LibraryConfig.build(
            library=lib, protocol=protocol, bam_file=bam, inputs=inputs_by_lib[lib], stranded_bam=stranded,
            disease_status=diseased, annotations=args.annotations, log=log,
            sample_size=args.genome_bins if protocol == PROTOCOL.GENOME else args.transcriptome_bins,
            distribution_fraction=args.distribution_fraction
        )
        libs.append(l)
    convert = {}
    for alias, command in args.external_conversion:
        if alias in convert:
            raise KeyError('duplicate alias names are not allowed', alias)
        convert[alias] = []
        open_option = False
        for item in re.split('\s+', command):
            if len(convert[alias]) > 0:
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
        stranded = str(TSV.tsv_boolean(stranded))
        SUPPORTED_TOOL.enforce(toolname)
        convert[alias] = ['convert_tool_output', inputfile, toolname, stranded]
    print(args)
    write_config(args.write, include_defaults=not args.no_defaults, libraries=libs, conversions=convert, log=log)


def time_diff(start, end):
    """
    >>> time_diff('2017-04-02 15:10:48.607195', '2017-04-03 17:00:32.671809')
    25.83
    >>> time_diff('2017-04-03 15:10:48.607195', '2017-04-03 17:00:32.671809')
    1.83
    """
    t1 = datetime.strptime(start, '%Y-%m-%d %H:%M:%S.%f')
    t2 = datetime.strptime(end, '%Y-%m-%d %H:%M:%S.%f')
    td = t1 - t2
    return td.seconds / 3600


def unique_exists(pattern, allow_none=False, get_newest=False):
    result = glob.glob(pattern)
    if len(result) == 1:
        return result[0]
    elif len(result) > 0:
        if get_newest:
            current_file = result[0]
            for f in result[1:]:
                stats1 = os.stat(current_file)
                stats2 = os.stat(f)
                if stats1.st_mtime < stats2.st_mtime:
                    current_file = f
            return current_file

        else:
            raise OSError('duplicate results:', result)
    elif allow_none:
        return None
    else:
        raise OSError('no result found', pattern)


def parse_log_details(log_file):
    details = None
    with open(log_file) as job_out:
        lines = job_out.readlines()
        if len(lines) == 0:
            details = argparse.Namespace(status='empty', message='empty log file', filename=log_file)
        elif 'error' in lines[-1].lower():
            details = argparse.Namespace(status='crash', message=lines[-1].strip(), filename=log_file)
        else:
            details = argparse.Namespace(
                status='incomplete', message=lines[-1].strip(), filename=log_file, last_modified=time.ctime(os.path.getmtime(log_file)))
    return details


def parse_runtime_from_log(log_file):

    with open(log_file) as job_out:
        lines = job_out.readlines()
        if len(lines) > 0:
            for line in lines[-10:]:
                m = re.match('^\s*run time \(s\): (\d+)\s*$', line)
                if m:
                    return int(m.group(1))
    raise OSError('error in parsing the log file for run times', log_file, '^\s*run time (s): (\d+)\s*$')


def convert_set_to_ranges(input_set):
    ranges = []
    for curr in sorted(list(input_set)):
        if ranges:
            if ranges[-1][1] + 1 == curr:
                ranges[-1] = (ranges[-1][0], curr)
                continue
        ranges.append((curr, curr))
    result = []
    for s, t in ranges:
        if s == t:
            result.append(str(s))
        else:
            result.append(str(s) + '-' + str(t))
    return ', '.join(result)


def check_library_dir(library_dir, verbose=False):
    """
    checks the completion of library directory and its clustering, validation and annotation subdirectories
    """
    log('Checking {} '.format(os.path.basename(library_dir)))

    clustering_dir = os.path.join(library_dir, 'clustering')
    if not os.path.isdir(clustering_dir):
        raise OSError('missing clustering directory', clustering_dir)

    clustering_complete_stamp = os.path.join(clustering_dir, 'MAVIS.COMPLETE')
    if not os.path.exists(clustering_complete_stamp):
        raise OSError('missing completion stamp', clustering_complete_stamp)
    log('clustering OK', time_stamp=False)
    # check for multiple batches as well as number of jobs run per batch
    batches = {}
    for job in [os.path.basename(p) for p in glob.glob(os.path.join(clustering_dir, 'batch*.tab'))]:
        batch = re.sub('-\d+.tab$', '', job)
        batches[batch] = batches.get(batch, 0) + 1
    if len(batches) > 1:
        raise OSError('Multiple clustering runs were found in the clustering directory', clustering_dir, list(batches.keys()))
    elif len(batches) == 0:
        raise OSError('No clustering output files were created', clustering_dir)

    batch_id, job_count = list(batches.items())[0]
    stamps = {'validation': {}, 'annotation': {}}
    run_times = {'validation': {}, 'annotation': {}}
    for stage_subdir in ['validation', 'annotation']:
        incomplete = set()
        missing_logs = set()
        crashes = {}
        curr_run_times = []
        printed_stage = False
        for job_task_id in [c for c in range(1, job_count + 1)]:
            stamp_pattern = os.path.join(library_dir, stage_subdir, batch_id + '-' + str(job_task_id), '*.COMPLETE')
            log_pattern = os.path.join(library_dir, stage_subdir, '*.' + str(job_task_id))
            try:
                stamp = unique_exists(stamp_pattern)
                stamps[stage_subdir][job_task_id] = os.path.getctime(stamp)
            except OSError:
                stamp = None
                if not printed_stage:
                    log(stage_subdir, 'FAIL', time_stamp=False)
                    printed_stage = True
                incomplete.add(job_task_id)
            try:
                logfile = unique_exists(log_pattern, get_newest=True)
            except OSError:
                logfile = None
                if not printed_stage:
                    log(stage_subdir, 'FAIL', time_stamp=False)
                    printed_stage = True
                missing_logs.add(job_task_id)

            if stamp and logfile:
                try:
                    rt = parse_runtime_from_log(logfile)
                except OSError as err:
                    incomplete.add(job_task_id)
                    crashes.setdefault(str(err), set()).add(job_task_id)
                else:
                    curr_run_times.append(rt)
                    run_times[stage_subdir][job_task_id] = rt
                    continue
            if logfile:
                if not printed_stage:
                    log(stage_subdir, 'FAIL', time_stamp=False)
                    printed_stage = True
                details = parse_log_details(logfile)
                if details.status == 'crash':
                    crashes.setdefault(details.message, set()).add(job_task_id)

        if incomplete or missing_logs or crashes:
            if incomplete:
                log(
                    '\t' + str(len(incomplete)), 'incomplete (jobs:', convert_set_to_ranges(incomplete) + ')', time_stamp=False)
            if missing_logs:
                log(
                    '\t' + str(len(missing_logs)), 'missing logs (jobs:', convert_set_to_ranges(missing_logs) + ')', time_stamp=False)
            if crashes:
                for msg, task_ids in crashes.items():
                    log('\t' + msg, '(jobs:', convert_set_to_ranges(task_ids) + ')', time_stamp=False)
        else:
            log(stage_subdir, 'OK', time_stamp=False)
            log('\trun time (s): {} (max), {} (total)'.format(max(curr_run_times), sum(curr_run_times)), time_stamp=False)

    max_rt = []
    for job, create_time in stamps['validation'].items():
        if create_time < os.path.getctime(clustering_complete_stamp):
            log('** ERROR: the clustering complete stamp must precede validation complete stamps')
        if job in stamps['annotation'] and stamps['annotation'][job] < create_time:
            log('** ERROR: A complete stamp for validation is after a complete stamp for annotation')
        if job in run_times['validation'] and job in run_times['annotation']:
            max_rt.append(run_times['validation'][job] + run_times['annotation'][job])

    stamp_times = list(stamps['validation'].values()) + list(stamps['annotation'].values())
    return (
        max(stamp_times) if stamp_times else None,
        max(max_rt) if len(max_rt) else None,
        sum(max_rt) if len(max_rt) else None
    )


def check_single_job(directory):
    name = os.path.basename(directory)
    stamp_pattern = os.path.join(directory, '*.COMPLETE')
    log_pattern = os.path.join(directory, '*.o*')
    logged_fail = False
    logfile = None
    stamp = None
    try:
        stamp = unique_exists(stamp_pattern)
    except OSError as err:
        log(name, 'FAIL')
        logged_fail = True
        log('\tINCOMPLETE: missing the complete stamp', stamp_pattern, time_stamp=False)
        return None, None
    try:
        logfile = unique_exists(log_pattern, allow_none=True, get_newest=True)
    except OSError as err:
        if not logged_fail:
            log(name, 'FAIL')
            logged_fail = True
        log('\tERROR:', err, time_stamp=False)

    rt = None
    if not stamp and not logfile:
        log('\t' + name, 'has not started', time_stamp=False)
    elif not stamp:
        log(parse_log_details(logfile), time_stamp=False)
    elif stamp:
        log(name, 'OK')
        rt = parse_runtime_from_log(logfile)
        log('\trun time (s):', rt, time_stamp=False)

    return (os.path.getctime(stamp) if stamp else None, rt)


def check_completion(target_dir):
    """
    Args:
        target_dir (str): path to the main pipeline output directory
    """
    if not os.path.isdir(target_dir):
        raise TypeError('expected a directory as input')
    library_dir_regex = '^[\w-]+_({})_({})$'.format('|'.join(DISEASE_STATUS.values()), '|'.join(PROTOCOL.values()))
    summary_stamp = None
    pairing_stamp = None
    library_stamps = []
    pipeline_total_rt = 0
    pipeline_max_rt = 0
    for subdir in sorted(glob.glob(os.path.join(target_dir, '*'))):
        name = os.path.basename(subdir)
        if name == 'summary':
            summary_stamp, rt = check_single_job(subdir)
            if rt:
                pipeline_total_rt += rt
                pipeline_max_rt += rt
        elif name == 'pairing':
            pairing_stamp, rt = check_single_job(subdir)
            if rt:
                pipeline_max_rt += rt
                pipeline_total_rt += rt
        elif re.match(library_dir_regex, name):
            last_timestamp, max_rt, total_rt = check_library_dir(subdir)
            if max_rt:
                pipeline_total_rt += total_rt
                pipeline_max_rt += max_rt
            if last_timestamp:
                library_stamps.append(last_timestamp)
        else:
            log('ignoring dir', subdir)
    library_stamp = max(library_stamps) if library_stamps else None

    if pairing_stamp and pairing_stamp < library_stamp:
        log('** ERROR: pairing completion stamp is older than library stamps')
    if summary_stamp and summary_stamp < pairing_stamp:
        log('** ERROR: summary stamp is older than pairing stamp')

    if pipeline_total_rt:
        log('pipeline total run time (s):', pipeline_total_rt)
    if pipeline_max_rt:
        log('pipeline max parallel run time (s):', pipeline_max_rt)


def main():
    def usage(err=None, detail=False):
        u = '\nusage: {} {{cluster,validate,annotate,pairing,summary,pipeline,config,checker}} [-h] [-v]'.format(PROGNAME)
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
        print(u)
        if detail:
            print(helpmenu)
        if err:
            print('{}: error:'.format(PROGNAME), err, '\n')
            exit(1)
        exit(0)

    start_time = int(time.time())

    if len(sys.argv) < 2:
        usage('the <pipeline step> argument is required')
    elif sys.argv[1] in ['-h', '--help']:
        usage(detail=True)
    elif sys.argv[1] in ['-v', '--version']:
        print('{} version {}'.format('mavis', __version__))
        exit(0)

    pstep = sys.argv.pop(1)
    sys.argv[0] = '{} {}'.format(sys.argv[0], pstep)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    augment_parser(required, optional, ['help', 'version'])

    if pstep == 'config':
        generate_config(parser, required, optional)
        exit(0)
    else:
        required.add_argument('-o', '--output', help='path to the output directory', required=True)
        if pstep == PIPELINE_STEP.PIPELINE:
            required.add_argument('config', help='path to the input pipeline configuration file')
            augment_parser(required, optional, [])
        elif pstep == PIPELINE_STEP.CLUSTER:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            augment_parser(
                required, optional,
                ['library', 'protocol', 'stranded_bam', 'disease_status'] +
                ['annotations', 'masking'] + [k for k in vars(CLUSTER_DEFAULTS)]
            )
        elif pstep == PIPELINE_STEP.VALIDATE:
            required.add_argument('-n', '--input', help='path to the input file', required=True)
            augment_parser(
                required, optional,
                ['library', 'protocol', 'bam_file', 'read_length', 'stdev_fragment_size', 'median_fragment_size'] +
                ['stranded_bam', 'annotations', 'reference_genome', 'aligner_reference', 'masking'] +
                [k for k in vars(VALIDATION_DEFAULTS)]
            )
        elif pstep == PIPELINE_STEP.ANNOTATE:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            augment_parser(
                required, optional,
                ['annotations', 'reference_genome', 'masking', 'max_proximity', 'template_metadata'] +
                [k for k in vars(ANNOTATION_DEFAULTS)] +
                [k for k in vars(ILLUSTRATION_DEFAULTS)]
            )
        elif pstep == PIPELINE_STEP.PAIR:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            optional.add_argument(
                '-f', '--product_sequence_files', nargs='+', help='paths to fasta files with product sequences',
                required=False, default=[])
            augment_parser(
                required, optional,
                ['annotations', 'max_proximity'] +
                [k for k in vars(PAIRING_DEFAULTS)]
            )
        elif pstep == PIPELINE_STEP.SUMMARY:
            required.add_argument('-n', '--inputs', nargs='+', help='path to the input files', required=True)
            augment_parser(
                required, optional,
                ['annotations', 'dgv_annotation', 'flanking_call_distance', 'split_call_distance',
                 'contig_call_distance', 'spanning_call_distance'] +
                [k for k in vars(SUMMARY_DEFAULTS)]
            )
        elif pstep == PIPELINE_STEP.CHECKER:

            args = parser.parse_args()
            check_completion(args.output)
            exit(0)
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
        print(repr(err))
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
    print('rargs', rargs, 'args', args)
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


if __name__ == '__main__':
    main()
