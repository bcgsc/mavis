import subprocess
import unittest
from unittest import mock

from mavis.schedule import scheduler as _scheduler
from mavis.schedule import constants as _constants
from mavis.schedule import job as _job
from mavis.constants import SUBCOMMAND


class TestParseQstat(unittest.TestCase):

    # TODO: single job running
    # TODO: batch job running
    # TODO: single job complete

    def test_single_job_complete(self):
        content = """
Job Id: 9.torque01.bcgsc.ca
    Job_Name = subtest.sh
    Job_Owner = creisle@torque01.bcgsc.ca
    resources_used.cput = 00:00:00
    resources_used.vmem = 346716kb
    resources_used.walltime = 00:01:00
    resources_used.mem = 3624kb
    resources_used.energy_used = 0
    job_state = C
    queue = batch
    server = torque01.bcgsc.ca
    Checkpoint = u
    ctime = Tue May 29 09:37:00 2018
    Error_Path = torque01.bcgsc.ca:/projects/trans_scratch/validations/workspa
        ce/creisle/temp/subtest.sh.e9
    exec_host = torque01.bcgsc.ca/0
    Hold_Types = n
    Join_Path = n
    Keep_Files = n
    Mail_Points = a
    mtime = Tue May 29 09:38:01 2018
    Output_Path = torque01.bcgsc.ca:/projects/trans_scratch/validations/worksp
        ace/creisle/temp/subtest.sh.o9
    Priority = 0
    qtime = Tue May 29 09:37:00 2018
    Rerunable = True
    Resource_List.walltime = 01:00:00
    Resource_List.nodes = 1
    Resource_List.nodect = 1
    session_id = 25438
    Variable_List = PBS_O_QUEUE=batch,PBS_O_HOME=/home/creisle,
        PBS_O_LOGNAME=creisle,
        PBS_O_PATH=/home/creisle/applications/node-v10.1.0-linux-x64/bin:/hom
        e/creisle/.npm-packages/bin:/home/creisle/bin:/home/creisle/applicatio
        ns/centos06/python-3.6.1/bin:/projects/tumour_char/analysis_scripts/bi
        n/pog:/gsc/software/linux-x86_64-centos6/git-2.12.0/bin/:/usr/local/bi
        n:/usr/local/sbin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/p
        rojects/trans_scratch/software/pipeline_commands/:/home/creisle/bin,
        PBS_O_MAIL=/var/spool/mail/creisle,PBS_O_SHELL=/bin/bash,
        PBS_O_LANG=en_US.UTF-8,
        PBS_O_WORKDIR=/projects/trans_scratch/validations/workspace/creisle/t
        emp,PBS_O_HOST=torque01.bcgsc.ca,PBS_O_SERVER=torque01.bcgsc.ca
    euser = creisle
    egroup = users
    queue_type = E
    comment = Job started on Tue May 29 at 09:37
    etime = Tue May 29 09:37:00 2018
    exit_status = 0
    submit_args = subtest.sh
    start_time = Tue May 29 09:37:01 2018
    start_count = 1
    fault_tolerant = False
    comp_time = Tue May 29 09:38:01 2018
    job_radix = 0
    total_runtime = 60.481239
    submit_host = torque01.bcgsc.ca
    init_work_dir = /projects/trans_scratch/validations/workspace/creisle/temp

    request_version = 1

        """
        rows = _scheduler.TorqueScheduler().parse_qstat(content)
        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual(_constants.JOB_STATUS.COMPLETED, row['status'])
        self.assertEqual('9.torque01.bcgsc.ca', row['job_ident'])
        self.assertEqual('subtest.sh', row['name'])
        self.assertIs(None, row['task_ident'])
        self.assertEqual('', row['status_comment'])

    def test_array_job(self):
        content = """
Job Id: 48[1].torque01.bcgsc.ca
    Job_Name = MA_mock-A47933_batch-JT3CUggKXNStHcoFXYaGR3-1
    Job_Owner = creisle@torque01.bcgsc.ca
    job_state = C
    queue = batch
    server = torque01.bcgsc.ca
    Checkpoint = u
    ctime = Tue May 29 18:27:33 2018
    depend = afterokarray:43[].torque01.bcgsc.ca
    Error_Path = torque01.bcgsc.ca:/projects/trans_scratch/validations/workspa
        ce/creisle/temp/test_submission/output_torque/mock-A47933_diseased_tra
        nscriptome/annotate/batch-JT3CUggKXNStHcoFXYaGR3-/job---.log-1
    Join_Path = oe
    Keep_Files = n
    Mail_Points = a
    mtime = Tue May 29 18:27:33 2018
    Output_Path = torque01.bcgsc.ca:/projects/trans_scratch/validations/worksp
        ace/creisle/temp/test_submission/output_torque/mock-A47933_diseased_tr
        anscriptome/annotate/batch-JT3CUggKXNStHcoFXYaGR3-/job---.log-1
    Priority = 0
    qtime = Tue May 29 18:27:33 2018
    Rerunable = True
    Resource_List.mem = 12000mb
    Resource_List.walltime = 16:00:00
    Resource_List.nodes = 1
    Resource_List.nodect = 1
    Variable_List = PBS_ARRAYID=1,PBS_O_QUEUE=batch,PBS_O_HOME=/home/creisle,
        PBS_O_LOGNAME=creisle,
        PBS_O_PATH=/home/creisle/git/mavis/venv/bin:/home/creisle/application
        s/node-v10.1.0-linux-x64/bin:/home/creisle/.npm-packages/bin:/home/cre
        isle/bin:/home/creisle/applications/centos06/python-3.6.1/bin:/project
        s/tumour_char/analysis_scripts/bin/pog:/gsc/software/linux-x86_64-cent
        os6/git-2.12.0/bin/:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr
        /bin:/usr/local/sbin:/usr/sbin:/projects/trans_scratch/software/pipeli
        ne_commands/:/home/creisle/bin,PBS_O_MAIL=/var/spool/mail/creisle,
        PBS_O_SHELL=/bin/bash,PBS_O_LANG=en_US.UTF-8,
        PBS_O_WORKDIR=/projects/trans_scratch/validations/workspace/creisle/t
        emp/test_submission,PBS_O_HOST=torque01.bcgsc.ca,
        PBS_O_SERVER=torque01.bcgsc.ca,
        MANPATH=/home/creisle/.npm-packages/share/man:/home/creisle/applicati
        ons/centos06/python-3.6.1/man:/usr/local/share/man:/usr/share/man/over
        rides:/usr/share/man,XDG_SESSION_ID=1340,HOSTNAME=torque01.bcgsc.ca,
        SHELL=/bin/bash,TERM=xterm-256color,HISTSIZE=1000,CLICOLOR=1,
        SSH_CLIENT=10.9.202.242 35994 22,TMPDIR=/var/tmp/,
        PYTHONUNBUFFERED=True,MAVIS_MIN_CLUSTERS_PER_FILE=2,
        NODE_OPTIONS=--trace-warnings,SSH_TTY=/dev/pts/0,USER=creisle,
        SVN_EDITOR=vim,LS_COLORS=di=34;01;47:mi=100;31;01:ln=36;01:ex=01;32,
        MAVIS_SCHEDULER=TORQUE,VIRTUAL_ENV=/home/creisle/git/mavis/venv,
        SACCT_FORMAT=jobid%-18\\,jobname%45\\,user%-8\\,reqmem\\,elapsed\\,state\\,
        MaxRSS\\,AveRSS\\,Partition,
        PATH=/home/creisle/git/mavis/venv/bin:/home/creisle/applications/node
        -v10.1.0-linux-x64/bin:/home/creisle/.npm-packages/bin:/home/creisle/b
        in:/home/creisle/applications/centos06/python-3.6.1/bin:/projects/tumo
        ur_char/analysis_scripts/bin/pog:/gsc/software/linux-x86_64-centos6/gi
        t-2.12.0/bin/:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/bin:/
        usr/local/sbin:/usr/sbin:/projects/trans_scratch/software/pipeline_com
        mands/:/home/creisle/bin,MAIL=/var/spool/mail/creisle,
        _=/usr/local/bin/qsub,
        PWD=/projects/trans_scratch/validations/workspace/creisle/temp/test_s
        ubmission,XMODIFIERS=@im=none,LANG=en_US.UTF-8,
        MODULEPATH=/usr/share/Modules/modulefiles:/etc/modulefiles,
        LOADEDMODULES=,
        NODE_PATH=/home/creisle/.npm-packages/lib/node_modules,
        SQUEUE_FORMAT=%.12i %9P %45j %.8u %.2t %.10M %.6D %.8m %.14l %.4c %.2
        0R %E,HISTCONTROL=ignoredups,MAVIS_MAX_FILES=1,HOME=/home/creisle,
        SHLVL=2,LOGNAME=creisle,
        PYTHONPATH=/home/creisle/applications/centos06/python-3.6.1/bin:,
        SSH_CONNECTION=10.9.202.242 35994 10.9.220.231 22,
        ORIENTDB_HOME=/home/creisle/applications/orientdb/orientdb-community-
        2.2.34,MODULESHOME=/usr/share/Modules,
        LESSOPEN=||/usr/bin/lesspipe.sh %s,BROWSER=/usr/bin/google-chrome,
        NPM_PACKAGES=/home/creisle/.npm-packages,
        XDG_RUNTIME_DIR=/run/user/1365,
        BASH_FUNC_module()=() {  eval `/usr/bin/modulecmd bash $*`\\
}
    euser = creisle
    egroup = users
    queue_type = E
    comment = Job 48[].torque01.bcgsc.ca deleted because its dependency of arr
        ay 43[].torque01.bcgsc.ca can never be satisfied
    etime = Tue May 29 18:27:33 2018
    exit_status = 271
    submit_args = -j oe -l mem=12000mb -l walltime=16:00:00 -V -W depend=after
        okarray:43[].torque01.bcgsc.ca -N MA_mock-A47933_batch-JT3CUggKXNStHco
        FXYaGR3 -o /projects/trans_scratch/validations/workspace/creisle/temp/
        test_submission/output_torque/mock-A47933_diseased_transcriptome/annot
        ate/batch-JT3CUggKXNStHcoFXYaGR3-/job---.log -t 1 /projects/trans_scra
        tch/validations/workspace/creisle/temp/test_submission/output_torque/m
        ock-A47933_diseased_transcriptome/annotate/submit.sh
    job_array_id = 1
    fault_tolerant = False
    job_radix = 0
    submit_host = torque01.bcgsc.ca
    init_work_dir = /projects/trans_scratch/validations/workspace/creisle/temp
        /test_submission
    request_version = 1

        """
        rows = _scheduler.TorqueScheduler().parse_qstat(content)
        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual('48[].torque01.bcgsc.ca', row['job_ident'])
        self.assertIs(1, row['task_ident'])

    # TODO: single job error
    # TODO: batch job error
    # TODO: single job exiting
    # TODO: batch job exiting


class TestCancel(unittest.TestCase):

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_single_job(self, patcher):
        sched = _scheduler.TorqueScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        patcher.assert_called_with(['qdel', '1234'])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_array_job(self, patcher):
        sched = _scheduler.TorqueScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        for task in job.task_list:
            self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
        patcher.assert_called_with(['qdel', '1234'])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_array_job_task(self, patcher):
        sched = _scheduler.TorqueScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job, task_ident='4')
        self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, job.status)
        for i, task in enumerate(job.task_list):
            if i == 3:
                self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
            else:
                self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, task.status)
        patcher.assert_called_with(['qdel', '1234', '-t', '4'])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_bad_command(self, patcher):
        patcher.side_effect = [subprocess.CalledProcessError(1, 'cmd')]
        sched = _scheduler.TorqueScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        sched.cancel(job)
        patcher.assert_called_with(['qdel', '1234'])
        self.assertNotEqual(_constants.JOB_STATUS.CANCELLED, job.status)


class TestSubmit(unittest.TestCase):

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_job(self, patcher):
        patcher.side_effect = ['141.torque01.bcgsc.ca\n']
        job = _job.Job(
            stage=SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            name='MV1',
            memory_limit=1,
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL,
            script='script.sh'
        )

        sched = _scheduler.TorqueScheduler()
        sched.submit(job)
        self.assertEqual('141.torque01.bcgsc.ca', job.job_ident)
        patcher.assert_called_with([
            'qsub', '-j', 'oe', '-q', 'all', '-l', 'mem=1mb',
            '-l', 'walltime=16:00:00', '-V', '-N', 'MV1',
            '-o', 'output_dir/job-$PBS_JOBNAME-$PBS_JOBID.log',
            '-m', 'abef', '-M', 'me@example.com', 'script.sh'
        ])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_job_with_job_deps(self, patcher):
        patcher.side_effect = ['141.torque01.bcgsc.ca\n']
        job = _job.Job(
            stage=SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            name='MV1',
            memory_limit=1,
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL,
            script='script.sh',
            dependencies=[
                _job.Job(
                    stage=SUBCOMMAND.VALIDATE,
                    output_dir='output_dir',
                    job_ident='1234.torque01.bcgsc.ca'
                ),
                _job.Job(
                    stage=SUBCOMMAND.VALIDATE,
                    output_dir='output_dir',
                    job_ident='54.torque01.bcgsc.ca'
                )]
        )

        sched = _scheduler.TorqueScheduler()
        sched.submit(job)
        self.assertEqual('141.torque01.bcgsc.ca', job.job_ident)
        patcher.assert_called_with([
            'qsub', '-j', 'oe', '-q', 'all', '-l', 'mem=1mb',
            '-l', 'walltime=16:00:00', '-V',
            '-W depend=afterok:1234.torque01.bcgsc.ca:54.torque01.bcgsc.ca',
            '-N', 'MV1',
            '-o', 'output_dir/job-$PBS_JOBNAME-$PBS_JOBID.log',
            '-m', 'abef', '-M', 'me@example.com', 'script.sh'
        ])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_job_with_mixed_deps(self, patcher):
        patcher.side_effect = ['141.torque01.bcgsc.ca\n']
        job = _job.Job(
            stage=SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            name='MV1',
            memory_limit=1,
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL,
            script='script.sh',
            dependencies=[
                _job.Job(
                    stage=SUBCOMMAND.VALIDATE,
                    output_dir='output_dir',
                    job_ident='1234.torque01.bcgsc.ca'
                ),
                _job.Job(
                    stage=SUBCOMMAND.VALIDATE,
                    output_dir='output_dir',
                    job_ident='54.torque01.bcgsc.ca'
                ),
                _job.TorqueArrayJob(
                    stage=SUBCOMMAND.VALIDATE,
                    output_dir='output_dir',
                    job_ident='99[].torque01.bcgsc.ca',
                    task_list=5
                )]
        )

        sched = _scheduler.TorqueScheduler()
        sched.submit(job)
        self.assertEqual('141.torque01.bcgsc.ca', job.job_ident)
        patcher.assert_called_with([
            'qsub', '-j', 'oe', '-q', 'all', '-l', 'mem=1mb',
            '-l', 'walltime=16:00:00', '-V',
            '-W depend=afterokarray:99[][5].torque01.bcgsc.ca,afterok:1234.torque01.bcgsc.ca:54.torque01.bcgsc.ca',
            '-N', 'MV1',
            '-o', 'output_dir/job-$PBS_JOBNAME-$PBS_JOBID.log',
            '-m', 'abef', '-M', 'me@example.com', 'script.sh'
        ])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_array(self, patcher):
        patcher.side_effect = ['142[].torque01.bcgsc.ca\n']
        job = _job.TorqueArrayJob(
            stage=SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            name='MV1',
            memory_limit=1,
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL,
            script='script.sh',
            task_list=[1, 2, 3, 6, 9]
        )

        sched = _scheduler.TorqueScheduler(concurrency_limit=2)
        sched.submit(job)
        self.assertEqual('142[].torque01.bcgsc.ca', job.job_ident)
        patcher.assert_called_with([
            'qsub', '-j', 'oe', '-q', 'all', '-l', 'mem=1mb',
            '-l', 'walltime=16:00:00', '-V', '-N', 'MV1',
            '-o', 'output_dir/job-$PBS_JOBNAME-$PBS_JOBID-$PBS_ARRAYID.log',
            '-m', 'abef', '-M', 'me@example.com',
            '-t', '1-3,6,9%2',
            'script.sh'
        ])
