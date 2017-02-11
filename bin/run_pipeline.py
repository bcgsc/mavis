from argparse import Namespace
from datetime import datetime
import errno
import os
from configparser import ConfigParser, ExtendedInterpolation
from structural_variant.validate import DEFAULTS as VDEFAULTS
import warnings
from structural_variant.constants import PROTOCOL
import argparse
from sv_merge import main as merge


QSUB_HEADER = """#!/bin/bash
#$ -V
#$ -N {name}
#$ -q {queue}
#$ -o {output}
#$ -l mem_free={memory}G,mem_token={memory}G,h_vmem={memory}G
#$ -j y"""

DEFAULTS = Namespace(
    min_clusters_per_file=50,
    max_files=100,
    cluster_clique_size=15,
    cluster_radius=20,
    **VDEFAULTS.__dict__,
    no_filter=False
)


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def mkdirp(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise


def read_config(filepath):
    """
    reads the configuration settings from the configuration file
    """
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read(filepath)

    LIBRARY_REQ_ATTR = ['protocol', 'bam_file', 'read_length', 'median_insert_size', 'stdev_isize', 'inputs']
    TYPE_CHECK = DEFAULTS.__dict__

    config = {
        'qsub': Namespace(memory=12, queue='transabyss.q'),
        'reference': Namespace(
            template_metadata=os.path.join(os.path.dirname(__file__), 'cytoBand.txt'),
            reference_genome='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
            annotations='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
            masking='/home/creisle/svn/svmerge/trunk/hg19_masked_regions.tsv'
        )
    }

    defaults = dict()
    defaults.update(DEFAULTS.__dict__)
    for attr, value in parser['DEFAULTS'].items():
        if attr == 'protocol':
            PROTOCOL.enforce(value)
        elif attr in TYPE_CHECK and type(TYPE_CHECK[attr]) != type(value):
            try:
                value = type(TYPE_CHECK[attr])(value)
                defaults[attr] = value
            except ValueError:
                defaults[attr] = value
                warnings.warn('type check failed for attr {} with value {}'.format(attr, repr(value)))
    
    library_sections = []

    for sec in parser.sections():
        section = dict()
        section.update(config.get(sec, Namespace()).__dict__)
        if sec == 'DEFAULTS':
            continue
        elif sec not in ['reference', 'qsub', 'visualization']:  # assume this is a library configuration
            library_sections.append(sec)
            for attr in LIBRARY_REQ_ATTR:
                if not parser.has_option(sec, attr):
                    raise KeyError(
                        'missing one or more required attribute(s) for the library section',
                        sec, attr, LIBRARY_REQ_ATTR)
            section.update(defaults)
            section['library'] = sec
    
        for attr, value in parser[sec].items():
            if attr == 'protocol':
                PROTOCOL.enforce(value)
            elif attr in TYPE_CHECK and type(TYPE_CHECK[attr]) != type(value):
                try:
                    value = type(TYPE_CHECK[attr])(value)
                except ValueError:
                    warnings.warn('type check failed for attr {} with value {}'.format(attr, repr(value)))
            elif attr in ['stdev_isize', 'median_insert_size', 'read_length']:
                try:
                    value = int(value)
                except ValueError:
                    value = float(value)
            elif attr == 'inputs':
                value = value.split(';') if value else []
            section[attr] = value

        n = Namespace(**section)
        config[sec] = n

    return config, library_sections


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='path to the configuration file')
    parser.add_argument(
        '--output', '-o', help='path to the output directory', default=os.path.join(os.getcwd(), 'output'))
    parser.add_argument(
        '-f', '--force_overwrite', help='will overwrite existing files/directories if they exist',
        default=False, action='store_true')
    args = parser.parse_args()

    if os.path.exists(args.output) and not args.force_overwrite:
        print('error: must specify the overwrite option or delete the existing output directory')
        parser.print_help()
        exit(1)
    config, library_sections = read_config(args.config)
    # read the config
    # set up the directory structure and run svmerge
    for library in library_sections:
        section = config[library]
        base = os.path.join(args.output, '{}_{}'.format(library, section.protocol))
        log('setting up the directory structure for', library, 'as', base)
        base = os.path.join(args.output, '{}_{}'.format(library, section.protocol))
        mkdirp(os.path.join(base, 'clustering'))
        mkdirp(os.path.join(base, 'validation'))
        mkdirp(os.path.join(base, 'annotation'))
        
        # run the merge
        log('clustering')
        section.__dict__.update(config['reference'].__dict__)
        merge(section)


    # run sv_merge
    # set up sv_validate jobs
    # set up sv_annotate jobs
    # set up sv_pair jobs
    """
        output_folder = '{}/validation/{}_{}'.format(args.output, lib, protocol)
        qsub_file = '{}/qsub.sh'.format(output_folder)
        with open(qsub_file, 'w') as fh:
            log('writing:', qsub_file)
            fh.write('#!/bin/sh\n')
            fh.write('#$ -t {}-{}\n'.format(1, len(jobs)))  # array job
            fh.write('#$ -V\n')  # copy environment variables
            fh.write('#$ -N svmV_{}\n'.format(lib[-5:]))
            fh.write('#$ -j y\n')
            fh.write('#$ -q {}\n'.format(args.queue))
            fh.write('#$ -o {}/log/\n'.format(output_folder))
            fh.write('#$ -l mem_free=12G,mem_token=12G,h_vmem=12G\n')
            fh.write('echo "Starting job: $SGE_TASK_ID"\n')
            fh.write('python sv_validate.py -n {}$SGE_TASK_ID.tab \\\n\t-o {} \\\n\t-b {} \\\n\t-l {} \\\n'.format(
                clusterset_file_prefix,
                output_folder,
                BAM_FILE_ARGS[lib],
                lib
            ))
            temp = EvidenceSettings()
            opt = []
            for param, val in args.__dict__.items():
                if hasattr(temp, param):
                    opt.append('--{} {}'.format(param, val))
            fh.write('\t{}\n'.format(' \\\n\t'.join(opt)))
            fh.write('echo "Job complete: $SGE_TASK_ID"\n'"""

if __name__ == '__main__':
    main()
