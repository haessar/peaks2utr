#!/usr/bin/env python3
import csv
from dataclasses import dataclass

import gffutils

# User parameter
maxDistanceFromTranscript = 200

gff_in = "Orig.PRFA01000011.gff"
gff_out = gff_in.split('.')[0] + '.db'

FORWARD_PEAKS_FILENAME = "forward_peaks.broadPeak"
REVERSE_PEAKS_FILENAME = "reverse_peaks.broadPeak"

db = gffutils.create_db(gff_in, gff_out, force=True)


@dataclass
class Peak:
    chr: str
    start: int
    end: int
    name: str
    score: int
    strand: str
    signalValue: float
    pValue: float
    qValue: float

    def __post_init__(self):
        # Enforce typing
        for k, v in self.__annotations__.items():
            self.__setattr__(k, v(self.__getattribute__(k)))


out = ''

strands = [(FORWARD_PEAKS_FILENAME, '+'), (REVERSE_PEAKS_FILENAME, '-')]

for peaks_filename, strand in strands:
    with open(peaks_filename, 'r') as f:
        peaks = csv.reader(f, delimiter="\t")
        for peak in peaks:
            peak = Peak(*peak)
            peak_range = range(peak.start, peak.end)
            features = list(db.region(
                seqid=peak.chr,
                start=peak.start - maxDistanceFromTranscript,
                end=peak.end + maxDistanceFromTranscript,
                strand=strand)
            )
            if features:
                if any([f for f in features if f.featuretype == 'three_prime_UTR']):  # Same for both strands
                    print("3' UTR already annotated for features near peak %s" % peak.name)
                    continue
                genes = [f for f in features if f.featuretype == 'gene']
                genes = list(reversed(genes)) if strand == '-' else genes
                for idx, gene in enumerate(genes):
                    gene_range = range(gene.start, gene.end)
                    if set(peak_range).issubset(gene_range):
                        print("Peak %s wholly contained within gene %s" % (peak.name, gene.id))
                        continue
                    if len(genes) > idx + 1:
                        next_gene = genes[idx + 1]
                        next_gene_range = range(next_gene.start, next_gene.end)
                        if set(peak_range).intersection(gene_range) and set(peak_range).intersection(next_gene_range):
                            print("Peak %s overlapping gene %s and gene %s" % (peak.name, gene.id, next_gene.id))
                            continue
                    if strand == "+" and peak.end > gene.end:
                        utr = {'start': gene.end, 'end': peak.end}
                    elif strand == "-" and peak.start < gene.start:
                        utr = {'start': peak.start, 'end': gene.start}
                    else:
                        utr = {}
                    if utr:
                        print("PEAK %s CORRESPONDS TO 3' UTR OF GENE %s" % (peak.name, gene.id))
                        print("-----> UTR = ({}, {})".format(utr['start'], utr['end']))
                        attrs = dict(gene.attributes)
                        attrs['ID'] = [gene.id + "_UTR"]
                        attrs['Parent'] = [gene.id]
                        attrs['colour'] = ['3']
                        f = gffutils.Feature(
                            seqid=gene.chrom,
                            source="3pUTR_annotation",
                            featuretype="three_prime_UTR",
                            start=utr['start'],
                            end=utr['end'],
                            score='.',
                            strand=gene.strand,
                            frame='.',
                            attributes=attrs
                        )
                        out += str(f) + '\n'

            else:
                print("No features found near peak %s" % peak.name)

with open('three_prime_UTRs.gff', 'w') as fout:
    fout.write(out)
