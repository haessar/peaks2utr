import argparse
import csv
import os
import os.path
import subprocess

import pysam
import MACS2
import gffutils

from . import constants, criteria, models
from .annotations import Annotations, annotate_utr_for_peak

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


def main():
    """
    The main function / pipeline for peaks2utr.
    """
    cache_dir = os.path.join(os.getcwd(), '.cache')
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    cached = lambda fn: os.path.join(cache_dir, fn)

    argparser = prepare_argparser()
    args = argparser.parse_args()

    gff_db = cached(os.path.splitext(args.GFF_IN)[0] + '.db')
    if not os.path.isfile(gff_db):
        gffutils.create_db(args.GFF_IN, gff_db, force=True)

    bam_basename = os.path.splitext(args.BAM_IN)[0]
    if not os.path.isfile(cached(bam_basename + '.forward.bam')):
        pysam.view("--threads", "6", "-b", "-F", "20", "-o", cached(bam_basename + '.forward.bam'), args.BAM_IN, catch_stdout=False)
    
    if not os.path.isfile(cached(bam_basename + '.reverse.bam')):
        pysam.view("--threads", "6", "-b", "-f", "16", "-o", cached(bam_basename + '.reverse.bam'), args.BAM_IN, catch_stdout=False)

    if not os.path.isfile(cached("forward_peaks.broadPeak")):
        p1 = subprocess.Popen(["macs2", "callpeak", "-t", cached(bam_basename + '.forward.bam'), "-n", "forward", "--nomodel", "--extsize", "100", "--broad", "--outdir", cache_dir])
    if not os.path.isfile(cached("reverse_peaks.broadPeak")):
        p2 = subprocess.Popen(["macs2", "callpeak", "-t", cached(bam_basename + '.reverse.bam'), "-n", "reverse", "--nomodel", "--extsize", "100", "--broad", "--outdir", cache_dir])
    try:
        exit_codes = [p.wait() for p in (p1, p2)]
    except NameError:
        pass

    db = gffutils.FeatureDB(gff_db)
    annotations = Annotations()
    total_peaks = 0
    for strand in ["forward", "reverse"]:
        with open(cached(strand + "_peaks.broadPeak"), 'r') as fin:
            peaks = list(csv.reader(fin, delimiter="\t"))
            total_peaks += len(peaks)
            for peak in peaks:
                peak = models.Peak(*peak)
                peak.strand = constants.STRAND_MAP.get(strand)
                annotate_utr_for_peak(db, annotations, peak, args.max_distance, args.override_utr, args.five_prime_ext)
    with open('three_prime_UTRs.gff', 'w') as fout:
        fout.writelines(annotations)
    with open('summary_stats.txt', 'w') as fstats:
        fstats.write(format_stats_line("Total peaks", total_peaks))
        fstats.write(format_stats_line("Total 3' UTRs annotated", len(annotations)))
        fstats.write(format_stats_line("Peaks with no nearby features", total_peaks, annotations.no_features_counter))
        fstats.write(format_stats_line("Peaks corresponding to an already annotated 3' UTR", total_peaks,
                                       len(criteria.assert_not_already_annotated.fails)))
        fstats.write(format_stats_line("Peaks contained within a feature", total_peaks,
                                       len(criteria.assert_not_a_subset.fails)))
        fstats.write(format_stats_line("Peaks corresponding to 5'-end of a feature", total_peaks,
                                       len(criteria.assert_3_prime_end_and_truncate.fails)))
