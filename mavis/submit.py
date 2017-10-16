from datetime import timedelta

from .constants import MavisNamespace
from .util import log, WeakMavisNamespace

STD_OPTIONS = ['memory_limit', 'queue', 'time_limit', 'import_env']

OPTIONS = WeakMavisNamespace()
OPTIONS.add('queue', '', cast_type=str, defn='the queue jobs are to be submitted to')
OPTIONS.add('memory_limit', 16000, defn='the maximum number of megabytes (MB) any given job is allowed')  # 16 GB
OPTIONS.add('import_env', True, defn='flag to import environment variables')
OPTIONS.add('time_limit', 10 * 60 * 60, defn='the time in seconds any given jobs is allowed')  # 10 hours
OPTIONS.add('validation_memory', 16000, defn='default memory limit (MB) for the validation stage')
OPTIONS.add('trans_validation_memory', 18000, defn='default memory limit (MB) for the validation stage (for transcriptomes)')
OPTIONS.add('annotation_memory', 12000, defn='default memory limit (MB) for the annotation stage')
OPTIONS.add('scheduler', 'SLURM', defn='The scheduler being used')


def build_dependency_string(command, delim, jobs):
    if isinstance(jobs, str):
        return command.format(jobs)
    return command.format(delim.join([str(j) for j in jobs]))


SCHEDULER = MavisNamespace(
    SGE=MavisNamespace(
        shebang='#!/bin/bash',
        submit='qsub -terse',
        option_prefix='#$',
        jobname='-N {}'.format,
        dependency=lambda x: build_dependency_string('-hold_jid {}', ',', x),
        queue='-q {}'.format,
        memory_limit=lambda x: '-l mem_free={0}G,mem_token={0}G,h_vmem={0}G'.format(x // 1000),
        join_output=lambda x: '-j {}'.format('y' if x else 'n'),
        import_env=lambda x: '-V',
        stdout='-o {}/sge-$JOB_NAME-$JOB_ID.log'.format,
        time_limit=lambda x: '-l h_rt={}'.format(str(timedelta(seconds=x)))
    ),
    SLURM=MavisNamespace(
        shebang='#!/bin/bash -l',
        submit='sbatch',
        option_prefix='#SBATCH',
        jobname='-J {}'.format,
        memory_limit='--mem {}M'.format,
        time_limit=lambda x: '-t {}'.format(str(timedelta(seconds=x))),
        stdout='-o {}/slurm-%x-%j.log'.format,
        dependency=lambda x: build_dependency_string('--dependency=afterok:{}', ':', x),
        import_env=lambda x: '--export=ALL',
        queue='--partition={}'.format
    )
)


class SubmissionScript:
    """
    holds scheduler options and build submissions scripts
    """
    def __init__(self, content, scheduler='SGE', **kwargs):
        self.options = {k: kwargs.pop(k, OPTIONS.get(k, None)) for k in ['jobname', 'stdout'] + STD_OPTIONS}
        if kwargs:
            raise TypeError('unexpected argument(s):', list(kwargs.keys()))
        self.scheduler = scheduler
        if scheduler not in SCHEDULER:
            raise ValueError('invalid scheduler', scheduler, 'expected', SCHEDULER.keys())
        for option, value in self.options.items():
            if value and option not in SCHEDULER[self.scheduler]:
                raise ValueError('scheduler', scheduler, 'does not support the option', option)
        self.content = content

    def __getattribute__(self, key):
        if key == 'options' or key not in self.options:
            return object.__getattribute__(self, key)
        return self.options[key]

    def build_header(self):
        """returns the header line detailing the scheduler-specific submission options"""
        config = SCHEDULER[self.scheduler]
        header = [config.shebang]
        if self.scheduler == 'SGE':
            header.append(config.option_prefix + ' ' + config['join_output'](True))
        for option, value in sorted(self.options.items()):
            if value is not None and value != '' and option in config:
                line = config.option_prefix + ' ' + config[option](value)
                header.append(line)
        return header

    def write(self, filepath):
        log('writing:', filepath)
        with open(filepath, 'w') as fh:
            for line in self.build_header():
                fh.write(line + '\n')
            fh.write('\n' + self.content + '\n')
        return filepath
