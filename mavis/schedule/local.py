import multiprocessing
import os

import shortuuid

from ..util import LOG

from .job import Job
from .scheduler import Scheduler
from .constants import JOB_STATUS


MAX_PROCESSES = 2

class LocalJob(Job):

    def __init__(self, args, func, rank=None, response=None, *pos, **kwargs):
        self.args = args
        self.func = func
        self.response = response
        self.rank = rank
        Job.__init__(self, *pos, **kwargs)

    def check_complete(self):
        return os.path.exists(self.complete_stamp())

    def flatten(self):
        result = Job.flatten(self)
        omit = {'script', 'rank', 'response', 'func', 'queue', 'import _env', 'mail_user', 'mail_type'}
        return {k:v for k, v in result.items() if k not in omit}


class LocalScheduler(Scheduler):
    """
    Scheduler class for dealing with running mavis locally
    """
    NAME = 'LOCAL'

    def __init__(self, concurrency_limit=None):
        """
        Args
            concurrency_limit (int): Size of the pool, the maximum allowed concurrent processes. Defaults to one less than the total number available
        """
        self.concurrency_limit = multiprocessing.cpu_count() - 1 if concurrency_limit is None else concurrency_limit
        self.pool = multiprocessing.Pool(self.concurrency_limit)
        self.submitted = {}  # submitted jobs process response objects by job ID

    def submit(self, job):
        """
        Add a job to the pool
        """
        if not self.pool:
            self.pool = multiprocessing.Pool(self.concurrency_limit)
        if not job.job_ident:
            job.job_ident = str(shortuuid.uuid())
            job.status = JOB_STATUS.SUBMITTED
        args = [arg.format(job_ident=job.job_ident, name=job.name) for arg in job.args]
        # if this job exists in the pool, return its response object
        if job.job_ident in self.submitted:
            return self.submitted[job.job_ident]
        # otherwise add it to the pool
        job.response = self.pool.apply_async(job.func, (args,))  # no arguments, defined all in the job object
        self.submitted[job.job_ident] = job
        job.rank = len(self.submitted)
        LOG('submitted', job.name, indent_level=1)
        return job

    def wait(self):
        """
        wait for everything in the current pool to finish
        """
        if not self.pool:
            return
        self.pool.close()
        self.pool.join()
        self.pool = None
        for job in self.submitted.values():
            self.update_info(job)

    def jobs_completed(self):
        return sum([1 for job in self.submitted.values() if job.status == JOB_STATUS.COMPLETED or job.response.ready()])

    def jobs_running(self):
        jobs = len(self.submitted) - self.jobs_completed()
        return min(self.concurrency_limit, jobs)

    def _check_running(self, job):
        return job.rank <= self.jobs_completed() + self.jobs_running()

    def update_info(self, job):
        """
        Args
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
            if job.response.ready():
                if job.response.successful():
                    job.status = JOB_STATUS.COMPLETED
                else:
                    job.status = JOB_STATUS.FAILED
            elif job.rank <= self.jobs_completed() + self.jobs_running():
                job.status = JOB_STATUS.RUNNING
            else:
                job.status = JOB_STATUS.PENDING
        else:
            job.status = JOB_STATUS.UNKNOWN
