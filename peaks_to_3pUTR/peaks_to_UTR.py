#!/usr/bin/env python3
from abc import ABC
import argparse
import collections
import csv

import gffutils

from peaks_to_3pUTR import criteria

FORWARD_PEAKS_FILENAME = "forward_peaks.broadPeak"
REVERSE_PEAKS_FILENAME = "reverse_peaks.broadPeak"

parser = argparse.ArgumentParser(
    description="Build an annotation GFF file containing three_prime_UTR features based on read coverage peaks."
)
parser.add_argument('gff_in', help="input 'canonical' GFF file")
parser.add_argument('forward_peaks', nargs='?', default=FORWARD_PEAKS_FILENAME,
                    help="forward strand peaks in BED 6+3 format")
parser.add_argument('reverse_peaks', nargs='?', default=REVERSE_PEAKS_FILENAME,
                    help="reverse strand peaks in BED 6+3 format")
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
        del attrs['ID']
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
        featuretype=['gene', 'transcript'])
    )
    genes = list(reversed(genes)) if peak.strand == '-' else genes
    if genes:
        for idx, gene in enumerate(genes):
            try:
                gene = genes[idx]
                gene.range = range(gene.start, gene.end)
                criteria.assert_not_already_annotated(db, peak, gene)
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


if __name__ == "__main__":
    args = parser.parse_args()

    gff_out = args.gff_in.split('.')[0] + '.db'
    db = gffutils.create_db(args.gff_in, gff_out, force=True)

    annotations = Annotations()
    strands = [(args.forward_peaks, '+'), (args.reverse_peaks, '-')]
    for peaks_filename, strand in strands:
        with open(peaks_filename, 'r') as fin:
            peaks = csv.reader(fin, delimiter="\t")
            for peak in peaks:
                peak = Peak(*peak)
                peak.strand = strand
                annotate_utr_for_peak(db, annotations, peak, args.max_distance)
    with open('three_prime_UTRs.gff', 'w') as fout:
        fout.writelines(annotations)
