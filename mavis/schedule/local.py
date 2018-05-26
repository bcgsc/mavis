import multiprocessing
import shlex
import shortuuid

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

    def __call__(self):
        self.job_ident = str(shortuuid.uuid())
        self.func(**self.args)


class LocalScheduler(Scheduler):
    NAME = 'LOCAL'

    def __init__(self, concurrency_limit=multiprocessing.cpu_count() - 1):
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
        job.response = self.pool.apply_async(job)  # no arguments, defined all in the job object
        self.submitted[job.job_ident] = job
        job.rank = len(self.submitted)

    def wait_on_jobs(self, job_list):
        if self.pool is None:
            self.pool = multiprocessing.Pool(self.concurrency_limit)
        for job in job_list:
            self.submit(job)
        self.pool.close()
        self.pool.join()
        self.pool = None
        for job in job_list:
            self.set_status(job)

    def submit_all(self, pipeline):
        # validations
        self.wait_on_jobs(pipeline.validations)
        # annotations
        self.wait_on_jobs(pipeline.annotations)
        # pairing
        self.wait_on_jobs(pipeline.pairing)
        # summary
        self.wait_on_jobs(pipeline.summary)

    def jobs_completed(self):
        return sum([1 for job in self.submitted.values() if job.response.ready()])

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
            if job.response.sucessful():
                job.status = JOB_STATUS.COMPLETE
            else:
                job.status = JOB_STATUS.FAILED
        elif job.rank <= self.jobs_completed() + self.jobs_running():
            job.status = JOB_STATUS.RUNNING
        else:
            job.status = JOB_STATUS.PENDING


def parse_script(filename):
    """
    Given some mavis bash script, parse and return just the mavis command
    """
    with open(filename, 'r') as fh:
        lines = [l.strip() for l in fh.readlines() if not l.startswith('#') and l]
        lines = ' '.join(lines)
        if sub:
            for original, replacement in sub:
                lines.replace(original, replacement)
        lex = shlex.split(lines)
        return [a.strip() for a in lex]
