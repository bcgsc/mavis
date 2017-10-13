from datetime import timedelta

from .constants import MavisNamespace
from .util import log, WeakMavisNamespace


OPTIONS = WeakMavisNamespace(
    jobname=None,
    dependency=None,
    queue=None,
    memory_limit=16000,  # 16 GB
    import_env=True,
    stdout=None,
    stderr=None,
    time_limit=20 * 60 * 60,  # 20 hours
    join_output=False
)


SCHEDULER = MavisNamespace(
    SGE=MavisNamespace(
        shebang='#!/bin/bash',
        submit='qsub',
        option_prefix='#$',
        jobname='-N {}'.format,
        dependency='-Wdepend=afterok:{}'.format,
        queue='-q {}'.format,
        memory_limit=lambda x: '-l mem_free={0}G,mem_token={0}G,h_vmem={0}G'.format(x // 1000),
        join_output=lambda x: '-j {}'.format('y' if x else 'n'),
        import_env=lambda x: '-V',
        stderr='-e {}'.format,
        stdout='-o {}'.format,
        time_limit=lambda x: '-l h_rt={}'.format(str(timedelta(seconds=x)))
    ),
    SLURM=MavisNamespace(
        shebang='#!/bin/bash -l',
        submit='sbatch',
        option_prefix='#SBATCH',
        jobname='-J {}'.format,
        memory_limit='--mem {}'.format,
        time_limit=lambda x: '-t {}'.format(str(timedelta(seconds=x))),
        stdout='-o {}'.format,
        stderr='-e {}'.format,
        dependency='--dependency=afterok:{}'.format,
        import_env=lambda x: '--export=ALL'
    )
)


class SubmissionScript:
    """
    holds scheduler options and build submissions scripts
    """
    def __init__(self, content, scheduler='SGE', **kwargs):
        self.options = {k: kwargs.pop(k, OPTIONS[k]) for k in OPTIONS}
        if not self.options['join_output']:
            self.options['join_output'] = None
        if not self.options['import_env']:
            self.options['import_env'] = None
        if kwargs:
            raise TypeError('unexpected argument(s):', list(kwargs.keys()))
        self.scheduler = scheduler
        if scheduler not in SCHEDULER:
            raise ValueError('invalid scheduler', scheduler, 'expected', SCHEDULER.keys())
        for option, value in self.options.items():
            if value is not None and value != OPTIONS[option] and option not in SCHEDULER[self.scheduler]:
                raise ValueError('scheduler', scheduler, 'does not support the option', option)
        if self.stderr and self.join_output:
            raise ValueError('stderr cannot be specified since join_output is set')
        self.content = content

    def __getattribute__(self, key):
        if key == 'options' or key not in self.options:
            return object.__getattribute__(self, key)
        return self.options[key]

    def build_header(self):
        """returns the header line detailing the scheduler-specific submission options"""
        config = SCHEDULER[self.scheduler]
        header = [config.shebang]
        for option, value in sorted(self.options.items()):
            if value is not None and option in config:
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
