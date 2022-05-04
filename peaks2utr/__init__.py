import argparse
import asyncio
import csv
import logging
import math
import multiprocessing
import os
import os.path
import shutil
import sys

import psutil
from tqdm import tqdm

from . import constants, models
from .annotations import Annotations, NoNearbyFeatures, batch_annotate_strand
from .utils import cached, iter_batches, yield_from_process, limit_memory
from .preprocess import BAMSplitter, call_peaks, create_db
from .postprocess import merge_and_gt_gff3_sort, write_summary_stats


def prepare_argparser():
    parser = argparse.ArgumentParser(
        description="""
        Use MACS2 to build forward and reverse peaks files for given .bam file.
        Iterate peaks through set of criteria to determine UTR viability, before annotating in .gff file.
        """
    )
    parser.add_argument('GFF_IN', help="input 'canonical' annotations file in gff or gtf format.")
    parser.add_argument('BAM_IN', help="input reads file in bam format.")
    parser.add_argument('--max-distance', type=int, default=200,
                        help='maximum distance in bases that UTR can be from a transcript.')
    parser.add_argument('--override-utr', action="store_true", help="ignore already annotated 3' UTRs in criteria.")
    parser.add_argument('--extend-utr', action="store_true", help="extend previously existing 3' UTR annotations where possible.")
    parser.add_argument('--five-prime-ext', type=int, default=0,
                        help='a peak within this many bases of a gene\'s 5\'-end should be assumed to belong to it.')
    parser.add_argument('--skip-soft-clip', action="store_true", help="skip the resource-intensive logic to pileup soft-clipped read edges.")
    parser.add_argument('--min-pileups', type=int, default=10, help='Minimum number of piled-up mapped reads for UTR cut-off.')
    parser.add_argument('--min-poly-tail', type=int, default=10, help='Minimum length of poly-A/T tail considered in soft-clipped reads.')
    parser.add_argument('-p', '--processors', type=int, default=1, help="How many processor cores to use.")
    parser.add_argument('-f', '-force', '--force', action="store_true", help="Overwrite outputs if they exist.")
    parser.add_argument('-o', '--output', help="output filename.")
    parser.add_argument('--keep-cache', action="store_true", help="Keep cached files on run completion.")
    return parser


def main():
    limit_memory(constants.PERC_ALLOCATED_VRAM * psutil.virtual_memory().total / 100)
    asyncio.run(_main())


async def _main():
    """
    The main function / pipeline for peaks2utr.
    """
    try:
        # Change root logger level from WARNING (default) to NOTSET in order for all messages to be delegated.
        logging.getLogger().setLevel(logging.NOTSET)

        # Add stdout handler, with level INFO.
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console.setFormatter(formatter)
        logging.getLogger().addHandler(console)

        if not os.path.exists(constants.LOG_DIR):
            logging.info("Make .log directory")
            os.mkdir(constants.LOG_DIR)

        # Add file handler, with level DEBUG.
        fileHandler = logging.FileHandler(filename=os.path.join(constants.LOG_DIR, '{}_debug.log'.format(__package__)), mode="w")
        fileHandler.setLevel(logging.DEBUG)
        fileHandler.setFormatter(formatter)
        logging.getLogger().addHandler(fileHandler)        

        if not os.path.exists(constants.CACHE_DIR):
            logging.info("Make .cache directory")
            os.mkdir(constants.CACHE_DIR)

        argparser = prepare_argparser()
        args = argparser.parse_args()

        gff_basename = os.path.basename(os.path.splitext(args.GFF_IN)[0])
        new_gff_fn = gff_basename + ".new.gff" if not args.output else args.output

        if os.path.exists(new_gff_fn) and not args.force:
            logging.error("%s already exists. Re-run with -f flag to force overwrite of output files. Aborting." % new_gff_fn)
            sys.exit(1)

        if args.override_utr and args.extend_utr:
            logging.error("Only one of --extend-utr and --override-utr can be used simultaneously. Aborting.")
            sys.exit(1)

        bam_basename = os.path.basename(os.path.splitext(args.BAM_IN)[0])
        bs = BAMSplitter(bam_basename, args)
        bs.split_strands()
        if not args.skip_soft_clip:
            bs.split_read_groups()
            bs.pileup_soft_clipped_reads()

        db, _, _ = await asyncio.gather(
            create_db(args.GFF_IN),
            call_peaks(bam_basename, "forward"),
            call_peaks(bam_basename, "reverse")
        )

        def parse_peaks(strand):
            with open(cached(strand + "_peaks.broadPeak"), 'r') as fin:
                peaks = [models.Peak(*peak) for peak in csv.reader(fin, delimiter="\t")]
                for peak in peaks:
                    peak.strand = constants.STRAND_MAP.get(strand)
                return peaks
                
        peaks = parse_peaks('forward') + parse_peaks('reverse')
        total_peaks = len(peaks)
        queue = multiprocessing.Queue()
        
        processes = [batch_annotate_strand(db, batch, queue, args) for batch in iter_batches(peaks, math.ceil(total_peaks/args.processors))]
        for p in processes:
            p.start()

        annotations = Annotations()
        with tqdm(total=total_peaks, desc=f'{"INFO": <8} Iterating over peaks to annotate 3\' UTRs.') as pbar:            
            for p in processes:
                for result in yield_from_process(queue, p, pbar):
                    if result:
                        if result is NoNearbyFeatures:
                            annotations.no_features_counter += 1
                        else:
                            annotations.update(result)

        merge_and_gt_gff3_sort(db, annotations, new_gff_fn, args)
        write_summary_stats(annotations, total_peaks)

        logging.info("%s finished successfully." % __package__)
        await asyncio.sleep(1)
        sys.exit(0)
    except KeyboardInterrupt:
        logging.error("User interrupted processing. Aborting.")
        sys.exit(130)
    finally:
        if not args.keep_cache:
            logging.info("Clearing cache.")
            shutil.rmtree(constants.CACHE_DIR)
