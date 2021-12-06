import argparse
import asyncio
import csv
import logging
import os
import os.path
import sys

from tqdm import tqdm

from . import constants, criteria, models
from .annotations import Annotations, annotate_utr_for_peak
from .utils import cached, format_stats_line, multiprocess_over_iterable
from .preprocess import call_peaks, create_db, pysam_strand_split


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
                        help='maximum distance in bases that UTR can be from a transcript')
    parser.add_argument('--override-utr', action="store_true", help="ignore already annotated 3' UTRs in criteria")
    parser.add_argument('--five-prime-ext', type=int, default=0,
                        help='a peak within this many bases of a gene\'s 5\'-end should be assumed to belong to it')
    return parser


async def main():
    """
    The main function / pipeline for peaks2utr.
    """
    # Change root logger level from WARNING (default) to NOTSET in order for all messages to be delegated.
    logging.getLogger().setLevel(logging.NOTSET)

    # Add stdout handler, with level INFO.
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

    if not os.path.exists(constants.LOG_DIR):
        logging.info("Make .log directory")
        os.mkdir(constants.LOG_DIR)

    # Add file handler, with level DEBUG.
    fileHandler = logging.FileHandler(filename=os.path.join(constants.LOG_DIR, 'debug.log'))
    fileHandler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fileHandler.setFormatter(formatter)
    logging.getLogger().addHandler(fileHandler)

    
    if not os.path.exists(constants.CACHE_DIR):
        logging.info("Make .cache directory")
        os.mkdir(constants.CACHE_DIR)

    argparser = prepare_argparser()
    args = argparser.parse_args()

    bam_basename = os.path.basename(os.path.splitext(args.BAM_IN)[0])
    multiprocess_over_iterable(['forward', 'reverse'], pysam_strand_split, [bam_basename, args])

    db, _, _ = await asyncio.gather(
        create_db(args.GFF_IN),
        call_peaks(bam_basename, "forward"),
        call_peaks(bam_basename, "reverse")
    )

    annotations = Annotations()
    total_peaks = 0
    for strand in ["forward", "reverse"]:
        with open(cached(strand + "_peaks.broadPeak"), 'r') as fin:
            peaks = list(csv.reader(fin, delimiter="\t"))
            total_peaks += len(peaks)
            logging.info("Iterating over %s strand peaks to annotate 3' UTRs." % strand)
            for peak in tqdm(peaks):
                peak = models.Peak(*peak)
                peak.strand = constants.STRAND_MAP.get(strand)
                annotate_utr_for_peak(db, annotations, peak, args.max_distance, args.override_utr, args.five_prime_ext)
    with open('three_prime_UTRs.gff', 'w') as fout:
        logging.info("Writing annotations to GFF output file.")
        fout.writelines(annotations)
    with open('summary_stats.txt', 'w') as fstats:
        logging.info("Writing summary statistics file.")
        fstats.write(format_stats_line("Total peaks", total_peaks))
        fstats.write(format_stats_line("Total 3' UTRs annotated", len(annotations)))
        fstats.write(format_stats_line("Peaks with no nearby features", total_peaks, annotations.no_features_counter))
        fstats.write(format_stats_line("Peaks corresponding to an already annotated 3' UTR", total_peaks,
                                       len(criteria.assert_not_already_annotated.fails)))
        fstats.write(format_stats_line("Peaks contained within a feature", total_peaks,
                                       len(criteria.assert_not_a_subset.fails)))
        fstats.write(format_stats_line("Peaks corresponding to 5'-end of a feature", total_peaks,
                                       len(criteria.assert_3_prime_end_and_truncate.fails)))
    logging.info("%s finished successfully." % __package__)