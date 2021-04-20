from snakemake.utils import validate
from snakemake.exceptions import WorkflowError
import os
from typing import List, Dict
import re
import json
import pandas as pd

CONTAINER = 'creisle/mavis:latest'

def output_dir(*paths):
    return os.path.join(config['output_dir'], *paths)

INITIALIZED_CONFIG = output_dir('config.json')


try:
    # TODO: replace with URL so that the user does not need a copy of the config schema
    validate(
        config,
        os.path.join(os.getcwd(), 'mavis/schemas/config.json')
    )
    for key in [
        "libraries",
        "reference.annotations",
        "output_dir"
    ]:
        if key not in config:
            raise ValueError(f'missing required property: {key}')
except Exception as err:
    short_msg = ' '.join(str(err).split('\n')[:2]) # these can get super long
    raise WorkflowError(short_msg)

libraries = sorted(list(config['libraries']))
VALIDATE_OUTPUT = output_dir('{library}/validate/batch-{job_id}/validation-passed.tab')
CLUSTER_OUTPUT = output_dir('{library}/cluster/batch-{job_id}.tab')

# create the cluster inputs and guess the cluster sizes
def count_total_rows(filenames):
    row_count = 0
    for filename in filenames:
        df = pd.read_csv(filename, sep='\t').drop_duplicates()
        row_count += df.shape[0]
    return row_count


for library in libraries:
    lib_config = config['libraries'][library]
    if 'total_batches' in lib_config:
        continue
    inputs = []
    for assignment in lib_config['assign']:
        if assignment in config['convert']:
            inputs.extend(config['convert'][assignment]['inputs'])
        else:
            inputs.append(assignment)

    # if not input by user, estimate the clusters based on the input files
    max_files = config['cluster.max_files']
    min_rows = config['cluster.min_clusters_per_file']
    total_rows = count_total_rows(inputs)

    if round(total_rows / max_files) >= min_rows:
        # use max number of jobs
        lib_config['total_batches'] = max_files
    else:
        lib_config['total_batches'] = total_rows // min_rows


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
    run:
        with open(output_dir('config.raw.json'), 'w') as fh:
            fh.write(json.dumps(config, sort_keys=True, indent='  '))


rule init_config:
    input: rules.copy_config.output
    output: INITIALIZED_CONFIG
    container: CONTAINER
    shell: 'mavis setup --config {input} --outputfile {output}'


rule convert:
    output: output_dir('converted_outputs/{alias}.tab')
    input: rules.init_config.output
    log: output_dir('converted_outputs/snakemake.{alias}.log.txt')
    params:
        file_type=lambda w: config['convert'][w.alias]['file_type'],
        strand_specific=lambda w: config['convert'][w.alias]['strand_specific'],
        assume_no_untemplated=lambda w: config['convert'][w.alias]['assume_no_untemplated'],
        input_files=lambda w: config['convert'][w.alias]['inputs']
    container: CONTAINER
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
    log: output_dir('snakemake.cluster.{library}.log.txt')
    container: CONTAINER
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
        log: output_dir('{library}/validate/snakemake.batch-{job_id}.log.txt')
        container: CONTAINER
        shell:
            'mavis validate --config {rules.init_config.output}'
                + ' --library {wildcards.library}'
                + ' --inputs {params.inputfile}'
                + ' --output {params.dirname}'
                + ' &> {log}'


rule annotate:
    input: rules.validate.output if not config['skip_stage.validate'] else rules.cluster.output
    output: stamp=output_dir('{library}/annotate/batch-{job_id}/MAVIS.COMPLETE'),
        result=output_dir('{library}/annotate/batch-{job_id}/annotations.tab')
    log: output_dir('{library}/annotate/snakemake.batch-{job_id}.log.txt')
    container: CONTAINER
    shell:
        'mavis annotate --config {rules.init_config.output}'
            + ' --library {wildcards.library}'
            + ' --inputs {input}'
            + ' --output ' + output_dir('{wildcards.library}/annotate/batch-{wildcards.job_id}')
            + ' &> {log}'


rule pairing:
    input: expand(rules.annotate.output.result, zip, library=libs_args, job_id=jobs_args)
    output: stamp=output_dir('pairing/MAVIS.COMPLETE'),
        result=output_dir('pairing/mavis_paired.tab')
    params:
        dirname=output_dir('pairing')
    log: output_dir('snakemake.pairing.log.txt')
    container: CONTAINER
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
    log: output_dir('snakemake.summary.log.txt')
    container: CONTAINER
    shell:
        'mavis summary --config {rules.init_config.output}'
            + ' --inputs {input}'
            + ' --output {params.dirname}'
            + ' &> {log}'
