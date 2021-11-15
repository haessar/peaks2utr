#!/usr/bin/env python3
import collections
import csv

import gffutils

# User parameter
maxDistanceFromTranscript = 200

gff_in = "Orig.PRFA01000011.gff"
gff_out = gff_in.split('.')[0] + '.db'

forward_broad_peaks = "forward_peaks.broadPeak"
reverse_broad_peaks = "reverse_peaks.broadPeak"

db = gffutils.create_db(gff_in, gff_out, force=True)

Peak = collections.namedtuple("Peak", 'chr start end name score strand signalValue pValue qValue')

out = ''

strands = [(forward_broad_peaks, '+'), (reverse_broad_peaks, '-')]

# TODO adapt for reverse strand as well
f = open(forward_broad_peaks, "r")
broad_peaks = csv.reader(f, delimiter="\t")
for peak in broad_peaks: 
    peak = Peak(*peak)
    features = list(db.region(
        seqid=peak.chr,
        start=int(peak.start) - maxDistanceFromTranscript,
        end=int(peak.end) + maxDistanceFromTranscript,
        strand='+')
    )
    if features:
        if any([f for f in features if f.featuretype == 'three_prime_UTR']):
            print("3' UTR already annotated for features near peak %s" % peak.name)
            continue
        genes = [f for f in features if f.featuretype == 'gene']
        for idx, gene in enumerate(genes):
            if gene.start < int(peak.start) and gene.end > int(peak.end):
                print("Peak %s wholly contained within gene %s" % (peak.name, gene.id))
                continue
            if len(genes) > idx + 1:
                if int(peak.start) < gene.end and int(peak.end) > genes[idx + 1].start:
                    print("Peak %s overlapping gene %s and gene %s" % (peak.name, gene.id, genes[idx + 1].id))
                    continue
            if int(peak.end) > gene.end:
                print("PEAK %s CORRESPONDS TO 3' UTR OF GENE %s" % (peak.name, gene.id))
                print("-----> UTR = (%s, %s)" % (gene.end, peak.end))
                attrs = dict(gene.attributes)
                attrs['ID'] = [gene.id + "_UTR"]
                attrs['Parent'] = [gene.id]
                attrs['colour'] = ['3']
                f = gffutils.Feature(
                    seqid=gene.chrom,
                    source="3pUTR_annotation",
                    featuretype="three_prime_UTR",
                    start=gene.end,
                    end=peak.end,
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
