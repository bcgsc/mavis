import multiprocessing
import shlex


MAX_PROCESSES = 2

class LocalJob(Job):
    def __init__(self, args, fn, *pos, **kwargs):
        self.args = args
        self.fn = fn
        Job.__init__(self, *pos, **kwargs)

    def __call__(self):
        self.fn(**self.args)


class LocalArrayJob(ArrayJob):
    def __init__(self, args, fn, *pos, **kwargs):
        self.args = args
        self.fn = fn
        ArrayJob.__init__(self, *pos, **kwargs)

    def __call__(self, task_ident):
        args = {}
        for arg, val in self.args.items():
            try:
                val = val.format(task_ident=task_ident)
            except AttributeError:
                pass
            args[arg] = val
        self.fn(**args)


class LocalScheduler(Scheduler):
    NAME = SCHEDULER.NONE

    def __init__(self, pipeline, concurrency_limit=2):
        self.pipeline = pipeline
        self.pool = multiprocessing.Pool(concurrency_limit)

    def submit()


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
