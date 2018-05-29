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


class LocalScheduler(Scheduler):
    """
    Scheduler class for dealing with running mavis locally
    """
    NAME = 'LOCAL'

    def __init__(self, concurrency_limit=multiprocessing.cpu_count() - 1):
        """
        Args
            concurrency_limit (int): Size of the pool, the maximum allowed concurrent processes. Defaults to one less than the total number available
        """
        self.concurrency_limit = concurrency_limit
        self.pool = multiprocessing.Pool(concurrency_limit)
        self.submitted = {}  # submitted jobs process response objects by job ID

    def submit(self, job):
        """
        Add a job to the pool
        """
        if not job.job_ident:
            job.job_ident = str(shortuuid.uuid())
        # if this job exists in the pool, return its response object
        if job.job_ident in self.submitted:
            return self.submitted[job.job_ident]
        # otherwise add it to the pool
        print(job.name, job.args, job.func)
        job.response = self.pool.apply_async(job.func, kwds=job.args)  # no arguments, defined all in the job object
        self.submitted[job.job_ident] = job
        job.rank = len(self.submitted)
        LOG('submitted', job.name, indent_level=1)

    def wait_on_jobs(self, job_list):
        if self.pool is None:
            self.pool = multiprocessing.Pool(self.concurrency_limit)
        for job in job_list:
            self.submit(job)
        self.pool.close()
        for job in job_list:
            job.response.get()
        self.pool = None
        for job in job_list:
            self.set_status(job)
            print('job status', job.name, job.status)

    def submit_all(self, pipeline):
        # validations
        LOG('submitting {} validate jobs'.format(len(pipeline.validations)), time_stamp=True)
        self.wait_on_jobs(pipeline.validations)
        # annotations
        LOG('submitting {} annotate jobs'.format(len(pipeline.annotations)), time_stamp=True)
        self.wait_on_jobs(pipeline.annotations)
        # pairing
        LOG('submitting the pairing job', time_stamp=True)
        self.wait_on_jobs([pipeline.pairing])
        # summary
        LOG('submitting the summary job', time_stamp=True)
        self.wait_on_jobs([pipeline.summary])

    def jobs_completed(self):
        return sum([1 for job in self.submitted.values() if job.status == JOB_STATUS.COMPLETED or job.response.ready()])

    def jobs_running(self):
        jobs = len(self.submitted) - self.jobs_completed()
        return min(self.concurrency_limit, jobs)

    def _check_running(self, job):
        return job.rank <= self.jobs_completed() + self.jobs_running()

    def set_status(self, job):
        if job.job_ident not in self.submitted:
            job.status = JOB_STATUS.NOT_SUBMITTED
            return
        if job.response.ready():  # is the job complete?
            if job.response.successful():
                job.status = JOB_STATUS.COMPLETE
            else:
                job.status = JOB_STATUS.FAILED
        elif job.rank <= self.jobs_completed() + self.jobs_running():
            job.status = JOB_STATUS.RUNNING
        else:
            job.status = JOB_STATUS.PENDING
