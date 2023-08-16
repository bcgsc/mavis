from snakemake.exceptions import WorkflowError
import os
import json
from mavis_config import (
    count_total_rows,
    get_library_inputs,
    get_singularity_bindings,
    guess_total_batches,
    validate_config,
)
from mavis_config.constants import SUBCOMMAND

# env variable mainly for CI/CD
CONTAINER = os.environ.get('SNAKEMAKE_CONTAINER', 'docker://bcgsc/mavis:v3.1.1')
MAX_TIME = 57600
DEFAULT_MEMORY_MB = 16000


if 'output_dir' not in config:
    raise WorkflowError('output_dir is a required property of the configfile')


def output_dir(*paths):
    return os.path.join(config['output_dir'], *paths)


INITIALIZED_CONFIG = output_dir('config.json')
LOG_DIR = output_dir('logs')

# external schedulers will not create the log dir if it does not already exist
if not os.path.exists(LOG_DIR):
    os.makedirs(LOG_DIR, exist_ok=True)

try:
    validate_config(config, stage=SUBCOMMAND.SETUP)
except Exception as err:
    short_msg = ' '.join(str(err).split('\n')[:2]) # these can get super long
    raise WorkflowError(short_msg)

# ADD bindings for singularity
workflow._singularity_args = f'-B {",".join(get_singularity_bindings(config))}'

libraries = sorted(list(config['libraries']))
VALIDATE_OUTPUT = output_dir('{library}/validate/batch-{job_id}/validation-passed.tab')
CLUSTER_OUTPUT = output_dir('{library}/cluster/batch-{job_id}.tab')


for library in libraries:
    if 'total_batches' in config['libraries'][library]:
        continue

    # if not input by user, estimate the clusters based on the input files
    config['libraries'][library]['total_batches'] = guess_total_batches(
        config, get_library_inputs(config, library)
    )


libs_args = []
jobs_args = []
for library in libraries:
    for job_id in range(1, config['libraries'][library]['total_batches'] + 1):
        libs_args.append(library)
        jobs_args.append(job_id)


rule all:
    input: output_dir('summary/MAVIS.COMPLETE')


rule copy_config:
    output: output_dir('config.raw.json')
    resources:
        time_limit=MAX_TIME,
        mem_mb=4000,
        cpus=1,
        log_dir=LOG_DIR
    log: os.path.join(LOG_DIR, 'copy_config.snakemake.log.txt')
    run:
        with open(output_dir('config.raw.json'), 'w') as fh:
            fh.write(json.dumps(config, sort_keys=True, indent='  '))

# adds the bam stats and default settings
rule init_config:
    input: rules.copy_config.output
    output: INITIALIZED_CONFIG
    container: CONTAINER
    log: os.path.join(LOG_DIR, 'init_config.snakemake.log.txt')
    resources:
        time_limit=MAX_TIME,
        mem_mb=DEFAULT_MEMORY_MB,
        cpus=1,
        log_dir=LOG_DIR
    shell: 'mavis setup --config {input} --outputfile {output}'


rule convert:
    output: output_dir('converted_outputs/{alias}.tab')
    input: rules.init_config.output
    log: os.path.join(LOG_DIR, 'convert.snakemake.{alias}.log.txt')
    params:
        file_type=lambda w: config['convert'][w.alias]['file_type'],
        strand_specific=lambda w: config['convert'][w.alias]['strand_specific'],
        assume_no_untemplated=lambda w: config['convert'][w.alias]['assume_no_untemplated'],
        input_files=lambda w: config['convert'][w.alias]['inputs']
    container: CONTAINER
    resources:
        time_limit=MAX_TIME,
        mem_mb=DEFAULT_MEMORY_MB,
        cpus=1,
        log_dir=LOG_DIR
    shell:
        'mavis convert --file_type {params.file_type}'
            + ' --strand_specific {params.strand_specific}'
            + ' --assume_no_untemplated {params.assume_no_untemplated}'
            + ' --inputs {params.input_files}'
            + ' --outputfile {output}'
            + ' &> {log}'


def get_cluster_inputs(w):
    conversions = config['convert']
    inputs = []
    for assignment in config['libraries'][w.library]['assign']:
        if assignment in conversions:
            inputs.extend(expand(rules.convert.output, alias=assignment))
        else:
            inputs.append(assignment)

    return inputs


rule cluster:
    input: files=get_cluster_inputs,
        config=rules.init_config.output
    output: directory(output_dir('{library}/cluster'))
    log: os.path.join(LOG_DIR, 'snakemake.cluster.{library}.log.txt')
    container: CONTAINER
    resources:
        time_limit=MAX_TIME,
        mem_mb=DEFAULT_MEMORY_MB,
        cpus=1,
        log_dir=LOG_DIR
    shell:
        'mavis cluster --config {input.config}'
            + ' --library {wildcards.library}'
            + ' --inputs {input.files}'
            + ' --output {output}'
            + ' &> {log}'


if not config['skip_stage.validate']:
    rule validate:
        input: rules.cluster.output
        params:
            dirname=lambda w: output_dir(f'{w.library}/validate/batch-{w.job_id}'),
            inputfile=lambda w: expand(CLUSTER_OUTPUT, library=[w.library], job_id=[w.job_id])
        output: VALIDATE_OUTPUT
        log: os.path.join(LOG_DIR, '{library}.validate.snakemake.batch-{job_id}.log.txt')
        container: CONTAINER
        resources:
            time_limit=MAX_TIME,
            mem_mb=18000,
            cpus=2,
            log_dir=LOG_DIR
        shell:
            'mavis validate --config {rules.init_config.output}'
                + ' --library {wildcards.library}'
                + ' --inputs {params.inputfile}'
                + ' --output {params.dirname}'
                + ' &> {log}'


def get_annotate_input_file(wildcards):
    if not config['skip_stage.validate']:
        return expand(rules.validate.output, library=[wildcards.library], job_id=[wildcards.job_id])
    return expand(CLUSTER_OUTPUT, library=[wildcards.library], job_id=[wildcards.job_id])


rule annotate:
    input: rules.validate.output if not config['skip_stage.validate'] else rules.cluster.output
    output: stamp=output_dir('{library}/annotate/batch-{job_id}/MAVIS.COMPLETE'),
        result=output_dir('{library}/annotate/batch-{job_id}/annotations.tab')
    params:
        inputfile=get_annotate_input_file
    log: os.path.join(LOG_DIR, '{library}.annotate.snakemake.batch-{job_id}.log.txt')
    container: CONTAINER
    resources:
        time_limit=MAX_TIME,
        mem_mb=DEFAULT_MEMORY_MB,
        cpus=2,
        log_dir=LOG_DIR
    shell:
        'mavis annotate --config {rules.init_config.output}'
            + ' --library {wildcards.library}'
            + ' --inputs {params.inputfile}'
            + ' --output ' + output_dir('{wildcards.library}/annotate/batch-{wildcards.job_id}')
            + ' &> {log}'


rule pairing:
    input: expand(rules.annotate.output.result, zip, library=libs_args, job_id=jobs_args)
    output: stamp=output_dir('pairing/MAVIS.COMPLETE'),
        result=output_dir('pairing/mavis_paired.tab')
    params:
        dirname=output_dir('pairing')
    log: os.path.join(LOG_DIR, output_dir('snakemake.pairing.log.txt'))
    container: CONTAINER
    resources:
        time_limit=MAX_TIME,
        mem_mb=DEFAULT_MEMORY_MB,
        cpus=1,
        log_dir=LOG_DIR
    shell:
        'mavis pairing --config {rules.init_config.output}'
            + ' --inputs {input}'
            + ' --output {params.dirname}'
            + ' &> {log}'


rule summary:
    input: rules.pairing.output.result,
    output: output_dir('summary/MAVIS.COMPLETE')
    params:
        dirname=output_dir('summary')
    log: os.path.join(LOG_DIR, 'snakemake.summary.log.txt')
    container: CONTAINER
    resources:
        time_limit=MAX_TIME,
        mem_mb=DEFAULT_MEMORY_MB,
        cpus=1,
        log_dir=LOG_DIR
    shell:
        'mavis summary --config {rules.init_config.output}'
            + ' --inputs {input}'
            + ' --output {params.dirname}'
            + ' &> {log}'
