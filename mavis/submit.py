from datetime import timedelta

from .constants import MavisNamespace
from .util import log, WeakMavisNamespace


SCHEDULER = MavisNamespace(SGE='SGE', SLURM='SLURM', __name__='~mavis.submit.SCHEDULER')
""":class:`~mavis.constants.MavisNamespace`: scheduler types

- :term:`SGE`
- :term:`SLURM`
"""

MAIL_TYPE = MavisNamespace(BEGIN='BEGIN', END='END', FAIL='FAIL', ALL='ALL', NONE='NONE')
"""
only supporting common mail type options between sge and slurm
https://slurm.schedmd.com/sbatch.html
http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
"""

STD_OPTIONS = ['memory_limit', 'queue', 'time_limit', 'import_env', 'mail_user', 'mail_type']

OPTIONS = WeakMavisNamespace(__name__='~mavis.submit.options')
""":class:`~mavis.constants.MavisNamespace`: submission options

- :term:`queue`
- :term:`memory_limit`
- :term:`import_env`
- :term:`time_limit`
- :term:`validation_memory`
- :term:`trans_validation_memory`
- :term:`annotation_memory`
- :term:`scheduler`

"""
OPTIONS.add('queue', '', cast_type=str, defn='the queue jobs are to be submitted to')
OPTIONS.add('memory_limit', 16000, defn='the maximum number of megabytes (MB) any given job is allowed')  # 16 GB
OPTIONS.add('import_env', True, defn='flag to import environment variables')
OPTIONS.add('time_limit', 16 * 60 * 60, defn='the time in seconds any given jobs is allowed')  # 16 hours
OPTIONS.add('validation_memory', 16000, defn='default memory limit (MB) for the validation stage')
OPTIONS.add('trans_validation_memory', 18000, defn='default memory limit (MB) for the validation stage (for transcriptomes)')
OPTIONS.add('annotation_memory', 12000, defn='default memory limit (MB) for the annotation stage')
OPTIONS.add('scheduler', SCHEDULER.SLURM, defn='The scheduler being used', cast_type=SCHEDULER)
OPTIONS.add('mail_type', '', cast_type=MAIL_TYPE.enforce, defn='When to notify the mail_user (if given)')
OPTIONS.add('mail_user', '', defn='User to send notifications to')


def _sge_mail_type(mail_type):
    """
    from: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
    `b'     Mail is sent at the beginning of the job.
    `e'     Mail is sent at the end of the job.
    `a'     Mail is sent when the job is aborted or rescheduled.
    `s'     Mail is sent when the job is suspended.
    `n'     No mail is sent.
    """
    if mail_type == MAIL_TYPE.NONE:
        return 'n'
    elif mail_type == MAIL_TYPE.BEGIN:
        return 'b'
    elif mail_type == MAIL_TYPE.END:
        return 'e'
    elif mail_type == MAIL_TYPE.FAIL:
        return 'as'
    else:
        return 'abes'


def _sge_mail_type_csv(mail_type_csv):
    mail_options = set()
    for opt in [o for o in mail_type_csv.split(',') if o]:
        MAIL_TYPE.enforce(opt)
        mail_options.update(_sge_mail_type(opt))
    if 'n' in mail_options and len(mail_options) > 1:
        raise ValueError('cannot specify NONE with other mail type options', mail_options)
    return ''.join(sorted(mail_options))


def _slurm_mail_type_csv(mail_type_csv):
    mail_options = set()
    for opt in [o for o in mail_type_csv.split(',') if o]:
        mail_options.add(MAIL_TYPE.enforce(opt))
    if MAIL_TYPE.NONE in mail_options and len(mail_options) > 1:
        raise ValueError('cannot specify NONE with other mail type options', mail_options)
    return ','.join(sorted(mail_options))


def build_dependency_string(command, delim, jobs):
    if isinstance(jobs, str):
        return command.format(jobs)
    return command.format(delim.join([str(j) for j in jobs]))


SCHEDULER_CONFIG = MavisNamespace(
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
        time_limit=lambda x: '-l h_rt={}'.format(str(timedelta(seconds=x))),
        mail_type=lambda x: '-m {}'.format(_sge_mail_type_csv(x)),
        mail_user='-M {}'.format
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
        queue='--partition={}'.format,
        mail_type=lambda x: '--mail-type={}'.format(_slurm_mail_type_csv(x)),
        mail_user='--mail-user={}'.format
    )
)


class SubmissionScript:
    """
    holds scheduler options and build submissions scripts
    """
    def __init__(self, content, scheduler, **kwargs):
        self.scheduler = scheduler
        self.options = {k: kwargs.pop(k, OPTIONS.get(k, None)) for k in ['jobname', 'stdout'] + STD_OPTIONS}
        if kwargs:
            raise TypeError('unexpected argument(s):', list(kwargs.keys()))
        if self.scheduler not in SCHEDULER_CONFIG:
            raise ValueError('invalid scheduler', self.scheduler, 'expected', SCHEDULER_CONFIG.keys())
        for option, value in self.options.items():
            if value and option not in SCHEDULER_CONFIG[self.scheduler]:
                raise ValueError('scheduler', self.scheduler, 'does not support the option', option)
        self.content = content

    def __getattribute__(self, key):
        if key == 'options' or key not in self.options:
            return object.__getattribute__(self, key)
        return self.options[key]

    def build_header(self):
        """returns the header line detailing the scheduler-specific submission options"""
        config = SCHEDULER_CONFIG[self.scheduler]
        header = [config.shebang]
        if self.scheduler == 'SGE':
            header.append(config.option_prefix + ' ' + config['join_output'](True))
        for option, value in sorted(self.options.items()):
            if option == 'mail_type' and not self.options.get('mail_user', None):  # ignore mail type if mail user option is not given
                continue
            if value is not None and value != '' and option in config:
                line = config.option_prefix + ' ' + config[option](value)
                header.append(line)
        return header

    def write(self, filepath):
        """
        write a submission script to the input path
        """
        log('writing:', filepath)
        with open(filepath, 'w') as fh:
            for line in self.build_header():
                fh.write(line + '\n')
            fh.write('\n' + self.content + '\n')
        return filepath
