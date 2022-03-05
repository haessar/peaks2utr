import multiprocessing
import os.path
from queue import Empty
import resource

from .constants import CACHE_DIR
from .exceptions import EXCEPTIONS_MAP


class Counter:
    seen = set()

    def __init__(self):
        self.val = multiprocessing.Value('i', 0)
        self.lock = multiprocessing.Lock()

    def add(self, peak):
        if peak not in self.seen:
            with self.lock:
                self.val.value += 1
                self.seen.add(peak)

    @property
    def value(self):
        with self.lock:
            return self.val.value


def cached(filename):
    return os.path.join(CACHE_DIR, filename)


async def consume_lines(pipe, log_file):
    """
    Asynchronously write lines in pipe to log file.
    """
    with open(log_file, 'bw') as f:
        while line := await pipe.readline():
            f.write(line)


def multiprocess_over_iterable(iterable, function, args):
    """
    Assign a multiprocessing Process to call function for every item in iterable, passing this item
    as the function's first argument.
    Start each process and wait for them all to finish before returning.
    """
    jobs = []
    for x in iterable:
        p = multiprocessing.Process(target=function, args=[x] + args)
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()
        if job.exitcode != 0:
            raise EXCEPTIONS_MAP.get(function.__name__, Exception)


def format_stats_line(msg, total, numerator=None):
    """
    Format given statistics message with optional percentage.
    """
    msg += ": "
    if numerator is None:
        msg += "{}\n".format(total)
    else:
        msg += "{} ({}%)\n".format(numerator, int(100 * numerator / total))
    return msg


def iter_batches(lst, n):
    """
    Yield successive n-sized chunks from lst.
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def yield_from_process(q, p, pbar=None):
    """
    Yield items in queue q while each process p is alive. This prevents program from locking up when queue
    gets too large.
    Pass an optional tqdm progress bar (pbar) to keep a single progress bar running over multiple processes.
    """
    while p.is_alive():
        p.join(timeout=1)
        while True:
            try:
                yield q.get(block=False)
                if pbar:
                    pbar.update()
            except Empty:
                break


def limit_memory(maxsize):
    """
    Limit total available memory globally to maxsize bytes. Will throw MemoryError if breached.
    """
    _, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard))


def filter_nested_dict(node, threshold):
    """
    For an n-nested dictionary, filter out integer leaves with a minimum threashold value.
    """
    if isinstance(node, int):
        if node >= threshold:
            return node
    else:
        dupe_node = {}
        for k, v in node.items():
            cur_node = filter_nested_dict(v, threshold)
            if cur_node:
                dupe_node[k] = cur_node
        return dupe_node or None
