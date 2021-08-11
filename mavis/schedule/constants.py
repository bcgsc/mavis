from ..constants import MavisNamespace
from ..util import WeakMavisNamespace


JOB_STATUS = MavisNamespace(
    "SUBMITTED",
    "COMPLETED",
    "ERROR",
    "RUNNING",
    "FAILED",
    "PENDING",
    "CANCELLED",
    NOT_SUBMITTED="NOT SUBMITTED",
    UNKNOWN="UNKNOWN",
    __name__="mavis.schedule.constants.JOB_STATUS",
)


def cumulative_job_state(states):
    """
    Given a set of states, return a single state based on the reporting priority
    """
    priority = [
        JOB_STATUS.ERROR,
        JOB_STATUS.FAILED,
        JOB_STATUS.CANCELLED,
        JOB_STATUS.NOT_SUBMITTED,
        JOB_STATUS.PENDING,
        JOB_STATUS.SUBMITTED,
        JOB_STATUS.RUNNING,
        JOB_STATUS.COMPLETED,
    ]
    for state in priority:
        if state in states:
            return state
    return JOB_STATUS.NOT_SUBMITTED


SCHEDULER = MavisNamespace(
    "SGE", "SLURM", "TORQUE", "LOCAL", __name__="mavis.schedule.constants.SCHEDULER"
)
"""MavisNamespace: scheduler types

- [LOCAL](/glossary/#LOCAL)
- [SGE](/glossary/#SGE)
- [SLURM](/glossary/#SLURM)
- [TORQUE](/glossary/#TORQUE)
"""

MAIL_TYPE = MavisNamespace(
    "BEGIN", "END", "FAIL", "ALL", "NONE", __name__="mavis.schedule.constants.MAIL_TYPE"
)
"""
When the scheduler should notify [mail_user](/configuration/settings/#mail_user) about a job

- ``ALL`` - All other options (except none)
- ``BEGIN`` - Send an email when the job starts
- ``END`` - Send an email when the job has terminated
- ``FAIL`` - Send an email if the job fails
- ``NONE`` - Do not send email
"""

STD_OPTIONS = [
    "memory_limit",
    "queue",
    "time_limit",
    "import_env",
    "mail_user",
    "mail_type",
    "priority",
]

OPTIONS = WeakMavisNamespace(__name__="mavis.schedule.constants.options")
"""MavisNamespace: submission options

- [annotation_memory](/configuration/settings/#annotation_memory)
- [concurrency_limit](/configuration/settings/#concurrency_limit)
- [import_env](/configuration/settings/#import_env)
- [mail_type](/configuration/settings/#mail_type)
- [mail_user](/configuration/settings/#mail_user)
- [memory_limit](/configuration/settings/#memory_limit)
- [queue](/configuration/settings/#queue)
- [remote_head_ssh](/configuration/settings/#remote_head_ssh)
- [scheduler](/configuration/settings/#scheduler)
- [time_limit](/configuration/settings/#time_limit)
- [trans_validation_memory](/configuration/settings/#trans_validation_memory)
- [validation_memory](/configuration/settings/#validation_memory)

"""
OPTIONS.add("annotation_memory", 12000, defn="default memory limit (MB) for the annotation stage")
OPTIONS.add("import_env", True, defn="flag to import environment variables")
OPTIONS.add(
    "mail_type", MAIL_TYPE.NONE, cast_type=MAIL_TYPE, defn="When to notify the mail_user (if given)"
)
OPTIONS.add("mail_user", "", defn="User(s) to send notifications to")
OPTIONS.add(
    "memory_limit", 16000, defn="the maximum number of megabytes (MB) any given job is allowed"
)  # 16 GB
OPTIONS.add("priority", "", defn="extra priority of the jobs to be submitted")
OPTIONS.add("queue", "", cast_type=str, defn="the queue jobs are to be submitted to")
OPTIONS.add("scheduler", SCHEDULER.SLURM, defn="The scheduler being used", cast_type=SCHEDULER)
OPTIONS.add(
    "time_limit", 16 * 60 * 60, defn="the time in seconds any given jobs is allowed"
)  # 16 hours
OPTIONS.add(
    "trans_validation_memory",
    18000,
    defn="default memory limit (MB) for the validation stage (for transcriptomes)",
)
OPTIONS.add("validation_memory", 16000, defn="default memory limit (MB) for the validation stage")
OPTIONS.add(
    "concurrency_limit",
    None,
    nullable=True,
    cast_type=int,
    defn="The concurrency limit for tasks in any given job array or the number of concurrent processes allowed for a local run",
)
OPTIONS.add("remote_head_ssh", "", cast_type=str, defn="ssh target for remote scheduler commands")
