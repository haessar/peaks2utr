from multiprocessing import Process
import os.path

from constants import CACHE_DIR


def cached(filename):
    return os.path.join(CACHE_DIR, filename)


async def consume_lines(pipe, log_file):
    with open(log_file, 'bw') as f:
        while line := await pipe.readline():
            f.write(line)


def multiprocess_over_iterable(iterable, function, args):
    jobs = []
    for x in iterable:
        p = Process(target=function, args=args + [x])
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()


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
