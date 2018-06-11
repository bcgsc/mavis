import atexit
from concurrent import futures
from datetime import datetime
import logging
import multiprocessing
import os

import shortuuid

from ..util import LOG
from ..annotate.file_io import REFERENCE_DEFAULTS, ReferenceFile

from .job import Job
from .scheduler import Scheduler
from .constants import JOB_STATUS, SCHEDULER


class LocalJob(Job):

    def __init__(self, args, func, rank=None, response=None, *pos, **kwargs):
        """
        Args:
            args (list): A list of arguments to passed to the function given
            func (callable): the function to be run
            rank (int): rank of the job within the pool
            response (:class:`~concurrent.futures.Future`): the result from the subprocess
        """
        self.args = args
        self.func = func
        self.response = response
        self.rank = rank
        for filetype in REFERENCE_DEFAULTS:
            setattr(self, filetype, kwargs.pop(filetype, None))
        Job.__init__(self, *pos, **kwargs)

    def check_complete(self):
        """
        check that the complete stamp associated with this job exists
        """
        return os.path.exists(self.complete_stamp())

    def flatten(self):
        result = Job.flatten(self)
        omit = {'script', 'rank', 'response', 'func', 'queue', 'import _env', 'mail_user', 'mail_type'}
        return {k: v for k, v in result.items() if k not in omit}


def write_stamp_callback(response):
    if response.exception() or response.cancelled() or response.running():
        return
    try:
        LOG('writing:', response.complete_stamp, time_stamp=True, indent_level=1)
        with open(response.complete_stamp, 'w') as fh:
            fh.write('end: {}\n'.format(int(datetime.timestamp(datetime.utcnow()))))
    except Exception as err:
        LOG('error writing the complete stamp', level=logging.CRITICAL, indent_level=1)
        raise err


class LocalScheduler(Scheduler):
    """
    Scheduler class for dealing with running mavis locally
    """
    NAME = SCHEDULER.LOCAL
    """:attr:`~mavis.schedule.constants.SCHEDULER`: the type of scheduler"""

    def __init__(self, *pos, **kwargs):
        Scheduler.__init__(self, *pos, **kwargs)
        self.concurrency_limit = multiprocessing.cpu_count() - 1 if not self.concurrency_limit else self.concurrency_limit
        self.pool = None  # set this at the first submission
        self.submitted = {}  # submitted jobs process response objects by job ID
        atexit.register(self.close)  # makes the pool 'auto close' on normal python exit

    def submit(self, job):
        """
        Add a job to the pool

        Args:
            job (LocalJob): the job to be submitted
        """
        if self.pool is None:
            self.pool = futures.ProcessPoolExecutor(max_workers=self.concurrency_limit)
        if not job.job_ident:
            job.job_ident = str(shortuuid.uuid())
            job.status = JOB_STATUS.SUBMITTED
        args = [arg.format(job_ident=job.job_ident, name=job.name) for arg in job.args]
        # if this job exists in the pool, return its response object
        if job.job_ident in self.submitted:
            return self.submitted[job.job_ident]

        # load any reference files not cached into the parent memory space
        for filetype in [f for f in REFERENCE_DEFAULTS.keys() if f != 'aligner_reference']:
            if getattr(job, filetype) is not None:
                ref = ReferenceFile(filetype, getattr(job, filetype))
                ref.load(verbose=False)
        # otherwise add it to the pool
        job.response = self.pool.submit(job.func, args)  # no arguments, defined all in the job object
        setattr(job.response, 'complete_stamp', job.complete_stamp())
        job.response.add_done_callback(write_stamp_callback)
        self.submitted[job.job_ident] = job
        job.rank = len(self.submitted)
        LOG('submitted', job.name, indent_level=1)
        return job

    def wait(self):
        """
        wait for everything in the current pool to finish
        """
        if self.pool is None:
            return
        self.pool.shutdown(True)
        self.pool = None
        for job in self.submitted.values():
            self.update_info(job)

    def update_info(self, job):
        """
        Args:
            job (LocalJob): the job to check and update the status for
        """
        # check if the job has been submitted already and completed or partially run
        if not job.job_ident:
            job.status = JOB_STATUS.NOT_SUBMITTED
        elif os.path.exists(job.complete_stamp()):
            job.status = JOB_STATUS.COMPLETED
        elif os.path.exists(job.logfile()) and job.job_ident not in self.submitted:
            job.status = JOB_STATUS.UNKNOWN
        elif job.job_ident in self.submitted:
            if job.response.done():
                excpt = job.response.exception()
                if excpt is None:
                    job.status = JOB_STATUS.COMPLETED
                else:
                    job.status = JOB_STATUS.FAILED
                    job.status_comment = str(excpt)
            elif job.response.running():
                job.status = JOB_STATUS.RUNNING
            else:
                job.status = JOB_STATUS.PENDING
        else:
            job.status = JOB_STATUS.UNKNOWN

    def close(self):
        if self.pool is not None:
            self.pool.shutdown()
            self.pool = None
