import argparse
import logging
import multiprocessing
import os.path
from queue import Empty
import re
import resource
import sqlite3

import pysam

from .constants import FeatureTypes, CACHE_DIR
from .exceptions import EXCEPTIONS_MAP
from .models import FeatureDB


class CustomArgumentParser(argparse.ArgumentParser):
    def parse_args(self, args=None, namespace=None):
        args = super().parse_args(args, namespace)
        if args.do_pseudo:
            self.add_pseudo_featuretypes()
        return args
    
    @staticmethod
    def add_pseudo_featuretypes():
        FeatureTypes.Gene.append("pseudogene")
        FeatureTypes.GffTranscript.append("pseudogenic_transcript")
        FeatureTypes.Exon.append("pseudogenic_exon")


class Falsey:

    def __bool__(self):
        return False


class Counter:
    seen = set()

    def __init__(self):
        self.val = multiprocessing.Value('i', 0)
        self.lock = multiprocessing.Lock()

    def __int__(self):
        return self.value

    def add(self, key):
        """
        Add key to global seen set. This Counter will only increment if key is not a duplicate in _any_ Counter.
        """
        if key not in self.seen:
            with self.lock:
                self.val.value += 1
                self.seen.add(key)

    @property
    def value(self):
        with self.lock:
            return self.val.value


def cached(filename):
    return os.path.join(CACHE_DIR, os.path.basename(filename))


def connect_db(db_path):
    db = sqlite3.connect(db_path, check_same_thread=False)
    return FeatureDB(db)


def index_bam_file(bam_file, processors):
    if not os.path.isfile(cached(bam_file + '.bai')):
        logging.info("Indexing %s." % bam_file)
        pysam.index("-@", str(processors), bam_file)


async def consume_lines(pipe, log_file):
    """
    Asynchronously write lines in pipe to log file.
    """
    with open(log_file, 'bw') as f:
        while line := await pipe.readline():
            f.write(line)


def multiprocess_over_dict(f, d):
    """
    Assign a multiprocessing Process to call function f for every key-value pair in d, passing this item
    as the function's first argument.
    Start each process and wait for them all to finish before returning.
    """
    jobs = []
    for input, output in d.items():
        p = multiprocessing.Process(target=f, args=(input, output))
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()
        if job.exitcode != 0:
            raise EXCEPTIONS_MAP.get(f.__name__, Exception)


def format_stats_line(msg, total, numerator=None):
    """
    Format given statistics message with optional percentage.
    """
    msg += ": "
    if numerator is None:
        msg += "{}\n".format(total)
    else:
        msg += "{} ({}%)\n".format(numerator, round(100 * numerator / total))
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
    resource.setrlimit(resource.RLIMIT_AS, (int(maxsize), hard))


def filter_nested_dict(node, threshold):
    """
    For an n-nested dictionary, filter out integer leaves with a minimum threshold value.
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


def sum_nested_dicts(d1, d2):
    """
    For an n-nested dictionary, sum numeric values in leaves with matching keys.
    """
    def sum(v1, v2):
        if v2 is None:
            return v1
        try:
            return v1 + v2
        except TypeError:
            return sum_nested_dicts(v1, v2)
    result = d2.copy()
    result.update({k: sum(v, d2.get(k)) for k, v in d1.items()})
    return result


def features_dict_for_gene(db, gene, transcript=None):
    """
    Return a dictionary containing gene and all its child features.
    Pass an optional transcript when 3' UTR has been annotated to allow extension.
    """
    features = {"gene": gene}
    for idx, f in enumerate(db.children(gene)):
        if f.id != gene.id:
            if transcript and f.id == transcript.id:
                features.update({"transcript": transcript})
            else:
                features.update({"feature_{}".format(idx): f})
    return features


def get_output_filename(args):
    gff_base, gff_ext = os.path.splitext(args.GFF_IN)
    gff_basename = os.path.basename(gff_base)
    args.gtf_in = True if "gtf" in gff_ext else False
    if not args.output:
        output_fn = gff_basename + ".new"
        output_fn += ".gtf" if args.gtf_out else ".gff3"
    else:
        if not args.gtf_out and args.output.endswith(".gtf"):
            args.gtf_out = True
        elif args.gtf_out and re.search(r".gff(3){0,1}$", args.output):
            logging.warning(f"""--gtf option has been submitted alongside {args.output} output filename.
                            This will lead to a GTF formatted file with a GFF extension. Please ensure
                            you consider the desired outcome.""")
        output_fn = args.output
    return output_fn
