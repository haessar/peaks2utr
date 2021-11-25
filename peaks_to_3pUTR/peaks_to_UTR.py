#!/usr/bin/env python3
from abc import ABC
import argparse
import collections
import csv

import gffutils

from peaks_to_3pUTR import criteria

strand_map = {
    'forward': '+',
    'reverse': '-',
}

parser = argparse.ArgumentParser(
    description="Build an annotation GFF file containing three_prime_UTR features based on read coverage peaks."
)
parser.add_argument('gff_db', help="gffutils.FeatureDB built from 'canonical' GFF file")
parser.add_argument('peaks', nargs='?', help="peaks in BED 6+3 format")
parser.add_argument('strand', choices=('forward', 'reverse'))
parser.add_argument('--max-distance', type=int, default=200,
                    help='maximum distance in bases that UTR can be from a transcript')


class RangeMixin(ABC):
    start: int
    end: int

    @property
    def range(self):
        return set(range(self.start, self.end))


class Peak(RangeMixin):
    """
    MACS2 peak in BED 6+3 format.
    """
    def __init__(self, *args):
        self.chr = str(args[0])
        self.start = int(args[1])
        self.end = int(args[2])
        self.name = str(args[3])
        self.score = int(args[4])
        self.strand = str(args[5])
        self.signalValue = float(args[6])
        self.pValue = float(args[7])
        self.qValue = float(args[8])

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, str(self.__dict__))


class UTR(RangeMixin):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.feature = None

    def __repr__(self):
        return "<%s: (%s, %s)>" % (self.__class__.__name__, self.start, self.end)

    def __eq__(self, other):
        return self.range == other.range

    def generate_feature(self, gene):
        attrs = dict(gene.attributes)
        attrs.pop('ID', None)
        attrs['Parent'] = [gene.id]
        attrs['colour'] = ['3']
        self.feature = gffutils.Feature(
            seqid=gene.chrom,
            source="3pUTR_annotation",
            featuretype="three_prime_UTR",
            start=self.start,
            end=self.end,
            score='.',
            strand=gene.strand,
            frame='.',
            attributes=attrs
        )

    def is_valid(self):
        return self.end > self.start


class Annotations(collections.UserDict):
    def __init__(self):
        super().__init__()
        self.no_features_counter = 0

    def __setitem__(self, gene, new_utr):
        existing_utr = self.get(gene)
        if existing_utr:
            if new_utr.range.issubset(existing_utr.range):
                return
        self.data[gene] = new_utr

    def __iter__(self):
        for _, utr in self.data.items():
            yield str(utr.feature) + '\n'


def annotate_utr_for_peak(db, annotations, peak, max_distance):
    genes = list(db.region(
        seqid=peak.chr,
        start=peak.start - max_distance,
        end=peak.end + max_distance,
        strand=peak.strand,
        featuretype=['gene', 'transcript', 'protein_coding_gene'])
    )
    genes = list(reversed(genes)) if peak.strand == '-' else genes
    if genes:
        for idx, gene in enumerate(genes):
            try:
                gene = genes[idx]
                gene.range = range(gene.start, gene.end)
                criteria.assert_not_already_annotated(peak, gene, db)
                criteria.assert_not_a_subset(peak, gene)
                utr = UTR(start=peak.start, end=peak.end)
                criteria.assert_3_prime_end_and_truncate(peak, gene, utr)
                if len(genes) > idx + 1:
                    next_gene = genes[idx + 1]
                    criteria.belongs_to_next_gene(peak, next_gene)
                    next_gene.range = range(next_gene.start, next_gene.end)
                    criteria.truncate_5_prime_end(peak, next_gene, utr)
            except criteria.CriteriaFailure as e:
                print("%s - %s" % (type(e).__name__, e))
            else:
                if utr.is_valid():
                    print("PEAK %s CORRESPONDS TO 3' UTR OF GENE %s" % (peak.name, gene.id))
                    print(utr)
                    utr.generate_feature(gene)
                    annotations[gene.id] = utr
    else:
        print("No features found near peak %s" % peak.name)
        annotations.no_features_counter += 1


def write_stats_line(file, msg, total, numerator=None):
    """
    Format given statistics message with optional percentage.
    """
    msg += ": "
    if numerator is None:
        msg += "{}\n".format(total)
    else:
        msg += "{} ({}%)\n".format(numerator, int(100 * numerator / total))
    file.write(msg)


if __name__ == "__main__":
    args = parser.parse_args()

    db = gffutils.FeatureDB(args.gff_db)
    annotations = Annotations()
    with open(args.peaks, 'r') as fin:
        peaks = list(csv.reader(fin, delimiter="\t"))
        total_peaks = len(peaks)
        for peak in peaks:
            peak = Peak(*peak)
            peak.strand = strand_map.get(args.strand)
            annotate_utr_for_peak(db, annotations, peak, args.max_distance)
    with open(args.strand + '_three_prime_UTRs.gff', 'w') as fout:
        fout.writelines(annotations)
    with open(args.strand + '_stats.txt', 'w') as fstats:
        write_stats_line(fstats, "Total peaks", total_peaks)
        write_stats_line(fstats, "Total 3' UTRs annotated", len(annotations))
        write_stats_line(fstats, "Peaks with no nearby features", total_peaks, annotations.no_features_counter)
        write_stats_line(fstats, "Peaks corresponding to an already annotated 3' UTR", total_peaks,
                         len(criteria.assert_not_already_annotated.fails))
        write_stats_line(fstats, "Peaks contained within a feature", total_peaks,
                         len(criteria.assert_not_a_subset.fails))
        write_stats_line(fstats, "Peaks corresponding to 5'-end of a feature", total_peaks,
                         len(criteria.assert_3_prime_end_and_truncate.fails))
