import os
import re
import sys


# local modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis.annotate import load_reference_genes, load_reference_genome, load_masking_regions, load_templates
from mavis import validate
from mavis import cluster
from mavis import pairing
from mavis import annotate

from mavis.validate.constants import VALIDATION_DEFAULTS
import mavis.pipeline.config as pconf
from mavis.pipeline.util import log, mkdirp
from mavis.constants import PROTOCOL

VALIDATION_PASS_SUFFIX = '.validation-passed.tab'

QSUB_HEADER = """#!/bin/bash
#$ -V
#$ -N {name}
#$ -q {queue}
#$ -o {output}
#$ -l mem_free={memory}G,mem_token={memory}G,h_vmem={memory}G
#$ -j y"""


def main_pipeline(args, configs):
    # read the config
    # set up the directory structure and run svmerge
    annotation_files = []
    annotation_jobs = []
    for sec in configs:
        base = os.path.join(args.output, '{}_{}'.format(sec.library, sec.protocol))
        log('setting up the directory structure for', sec.library, 'as', base)
        base = os.path.join(args.output, '{}_{}'.format(sec.library, sec.protocol))
        cluster_output = mkdirp(os.path.join(base, 'clustering'))
        validation_output = mkdirp(os.path.join(base, 'validation'))
        annotation_output = mkdirp(os.path.join(base, 'annotation'))

        # run the merge
        log('clustering')
        merge_args = {}
        merge_args.update(sec.__dict__)
        merge_args.update(args.__dict__)
        merge_args['output'] = os.path.join(base, 'clustering')
        output_files = cluster.main(**merge_args)
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
            'output': validation_output,
            'masking': args.masking_filename,
            'reference_genome': args.reference_genome_filename,
            'blat_2bit_reference': args.blat_2bit_reference,
            'annotations': args.annotations_filename,
            'library': sec.library,
            'bam_file': sec.bam_file,
            'protocol': sec.protocol,
            'read_length': sec.read_length,
            'stdev_fragment_size': sec.stdev_fragment_size,
            'median_fragment_size': sec.median_fragment_size,
            'stranded_bam': sec.stranded_bam,
            'protocol': sec.protocol
        }
        for attr in sorted(VALIDATION_DEFAULTS.__dict__.keys()):
            validation_args[attr] = getattr(sec, attr)

        qsub = os.path.join(validation_output, 'qsub.sh')
        validation_jobname = 'validation_{}_{}'.format(sec.library, sec.protocol)
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=args.queue, memory=args.validate_memory, name=validation_jobname, output=validation_output
                ) + '\n')
            fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            temp = [
                '--{} {}'.format(k, v) for k, v in validation_args.items() if not isinstance(v, str) and v is not None]
            temp.extend(
                ['--{} "{}"'.format(k, v) for k, v in validation_args.items() if isinstance(v, str) and v is not None])
            validation_args = temp
            validation_args.append('-n {}$SGE_TASK_ID.tab'.format(merge_file_prefix))
            fh.write('python {} validate {}\n'.format(__file__, ' \\\n\t'.join(validation_args)))

        # set up the annotations job
        # for all files with the right suffix
        annotation_args = {
            'output': annotation_output,
            'reference_genome': args.reference_genome_filename,
            'annotations': args.annotations_filename,
            'template_metadata': args.template_metadata_filename,
            'min_orf_size': sec.min_orf_size,
            'max_orf_cap': sec.max_orf_cap,
            'min_domain_mapping_match': sec.min_domain_mapping_match,
            'domain_regex_filter': sec.domain_regex_filter,
            'max_proximity': sec.max_proximity
        }
        temp = [
            '--{} {}'.format(k, v) for k, v in annotation_args.items() if not isinstance(v, str) and v is not None]
        temp.extend(
            ['--{} "{}"'.format(k, v) for k, v in annotation_args.items() if isinstance(v, str) and v is not None])
        annotation_args = temp
        annotation_args.append('--inputs {}/*{}*{}'.format(
            validation_output, os.path.basename(merge_file_prefix), VALIDATION_PASS_SUFFIX))
        qsub = os.path.join(annotation_output, 'qsub.sh')
        annotation_jobname = 'annotation_{}_{}'.format(sec.library, sec.protocol)
        annotation_jobs.append(annotation_jobname)
        annotation_files.append(os.path.join(annotation_output, 'annotations.tab'))
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=args.queue, memory=args.default_memory, name=annotation_jobname, output=annotation_output
                ) + '\n')
            fh.write('#$ -hold_jid {}\n'.format(validation_jobname))
            fh.write('python {} annotate {}\n'.format(__file__, ' \\\n\t'.join(annotation_args)))

    # set up scripts for the pairing held on all of the annotation jobs
    pairing_output = mkdirp(os.path.join(args.output, 'pairing'))
    pairing_args = dict(
        output=pairing_output,
        split_call_distance=args.split_call_distance,
        contig_call_distance=args.contig_call_distance,
        flanking_call_distance=args.flanking_call_distance,
        max_proximity=args.max_proximity,
        annotations=args.annotations_filename,
        low_memory=args.low_memory
    )
    temp = ['--{} {}'.format(k, v) for k, v in pairing_args.items() if not isinstance(v, str) and v is not None]
    temp.extend(['--{} "{}"'.format(k, v) for k, v in pairing_args.items() if isinstance(v, str) and v is not None])
    temp.append('--inputs {}'.format(' '.join(annotation_files)))
    pairing_args = temp
    qsub = os.path.join(pairing_output, 'qsub.sh')
    with open(qsub, 'w') as fh:
        log('writing:', qsub)
        fh.write(
            QSUB_HEADER.format(
                queue=args.queue, memory=args.default_memory, name='mavis_pairing', output=pairing_output
            ) + '\n')
        fh.write('#$ -hold_jid {}\n'.format(' '.join(annotation_jobs)))
        fh.write('python {} pairing {}\n'.format(__file__, ' \\\n\t'.join(pairing_args)))


def main():
    def usage(err):
        name = os.path.basename(__file__)
        print('usage: {} {{cluster,validate,annotate,pairing,summary,pipeline}}'.format(name))
        print('{}: error:'.format(name), err)
        exit(1)

    if len(sys.argv) < 2:
        usage('the <pipeline step> argument is required')

    pstep = sys.argv.pop(1)
    sys.argv[0] = '{} {}'.format(sys.argv[0], pstep)

    if pstep == pconf.PIPELINE_STEP.SUMMARY:
        raise NotImplementedError('summary script has not been written')
    args = pconf.parse_arguments(pstep)
    config = []
    log('input arguments')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', repr(val), time_stamp=False)
    if pstep == pconf.PIPELINE_STEP.PIPELINE:
        if args.write:
            log('writing:', args.config)
            pconf.write_config(args.config, include_defaults=True)
            exit()
        else:
            temp, config = pconf.read_config(args.config)
            args.__dict__.update(temp.__dict__)
            for sec in config:
                sec.output = args.output
    # load the reference files if they have been given and reset the arguments to hold the original file name and the
    # loaded data
    if any([
        pstep not in [pconf.PIPELINE_STEP.PIPELINE, pconf.PIPELINE_STEP.VALIDATE],
        hasattr(args, 'uninformative_filter') and args.uninformative_filter,
        pstep in [pconf.PIPELINE_STEP.VALIDATE, pconf.PIPELINE_STEP.CLUSTER] and args.protocol == PROTOCOL.TRANS,
        pstep == pconf.PIPELINE_STEP.PIPELINE and any([sec.protocol == PROTOCOL.TRANS for sec in config])
    ]):
        log('loading:', args.annotations)
        args.annotations_filename = args.annotations
        args.annotations = load_reference_genes(args.annotations)
    else:
        args.annotations_filename = args.annotations
        args.annotations = None

    try:
        if pstep != pconf.PIPELINE_STEP.PIPELINE:
            log('loading:' if not args.low_memory else 'indexing:', args.reference_genome)
            args.reference_genome_filename = args.reference_genome
            args.reference_genome = load_reference_genome(args.reference_genome, args.low_memory)
        else:
            args.reference_genome_filename = args.reference_genome
            args.reference_genome = None
    except AttributeError:
        pass
    try:
        log('loading:', args.masking)
        args.masking_filename = args.masking
        args.masking = load_masking_regions(args.masking)
    except AttributeError as err:
        pass
    try:
        if pstep != pconf.PIPELINE_STEP.PIPELINE:
            log('loading:', args.template_metadata)
            args.template_metadata_filename = args.template_metadata
            args.template_metadata = load_templates(args.template_metadata)
        else:
            args.template_metadata_filename = args.template_metadata
            args.template_metadata = None
    except AttributeError:
        pass

    # decide which main function to execute
    if pstep == pconf.PIPELINE_STEP.CLUSTER:
        cluster.main(**args.__dict__)
    elif pstep == pconf.PIPELINE_STEP.VALIDATE:
        validate.main(**args.__dict__)
    elif pstep == pconf.PIPELINE_STEP.ANNOTATE:
        annotate.main(**args.__dict__)
    elif pstep == pconf.PIPELINE_STEP.PAIR:
        pairing.main(**args.__dict__)
    elif pstep == pconf.PIPELINE_STEP.SUMMARY:
        pass    # main_summary(args)
    else:  # PIPELINE
        main_pipeline(args, config)

if __name__ == '__main__':
    main()
