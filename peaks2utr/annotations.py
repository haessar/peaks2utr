import collections
import copy
import json
import logging
import multiprocessing
import sqlite3

import gffutils

from . import constants, criteria, models
from .utils import cached


class NoNearbyFeatures:
    pass


class Annotations(collections.UserDict):
    def __init__(self):
        super().__init__()
        self.no_features_counter = 0

    def __setitem__(self, gene, new_features):
        existing_features = self.get(gene)
        if existing_features:
            if new_features["utr"].range.issubset(existing_features["utr"].range):
                return
        self.data[gene] = new_features

    def __iter__(self):
        for _, features in self.data.items():            
            nf = features.copy()
            if "utr" in nf:
                nf["utr"] = nf["utr"].feature
            yield '\n'.join([str(f) for _, f in nf.items()]) + '\n'


def set_gene_range(gene, strand, five_prime_ext=0):
    """
    Assign a range to gene taking into account assumed 5' extension.
    """
    if strand == "+":
        gene.start -= five_prime_ext
    else:
        gene.end += five_prime_ext
    gene.range = range(gene.start, gene.end)


def _iter_peaks(db, peaks_batch, queue, args):
    for peak in peaks_batch:
        annotate_utr_for_peak(db, queue, peak, args.max_distance, args.override_utr, args.extend_utr, args.five_prime_ext)         


def batch_annotate_strand(db, peaks_batch, queue, args):
    """
    Create multiprocessing Process to handle batch of peaks. Connect to sqlite3 db for each batch to prevent
    serialization issues.
    """
    db = sqlite3.connect(db, check_same_thread=False)
    db = gffutils.FeatureDB(db)
    return multiprocessing.Process(target=_iter_peaks, args=(db, peaks_batch, queue, args))


def annotate_utr_for_peak(db, queue, peak, max_distance, override_utr=False, extend_utr=False, five_prime_ext=0):
    """
    Find genes in region of given peak and apply criteria to determine if 3' UTR exists for each.
    If so, add to multiprocessing Queue.
    """
    utr_found = False
    genes = list(db.region(
        seqid=peak.chr,
        start=peak.start - max_distance,
        end=peak.end + max_distance,
        strand=peak.strand,
        featuretype=['gene', 'transcript', 'protein_coding_gene'])
    )
    for k, v in constants.STRAND_MAP.items():
        if peak.strand == v:
            with open(cached(k + "_unmapped.json"), "r") as f:
                spat_pileups = json.load(f) or {}
    genes = sorted(genes, key=lambda x: x.start, reverse=False if peak.strand == "+" else True)
    if genes:
        for idx, gene in enumerate(genes):
            try:
                gene = copy.deepcopy(genes[idx])
                set_gene_range(gene, peak.strand)
                criteria.assert_whether_utr_already_annotated(peak, gene, db, override_utr, extend_utr)
                criteria.assert_not_a_subset(peak, gene)
                utr = models.UTR(start=peak.start, end=peak.end)
                criteria.assert_3_prime_end_and_truncate(peak, gene, utr)
                if len(genes) > idx + 1:
                    next_gene = copy.deepcopy(genes[idx + 1])
                    set_gene_range(next_gene, peak.strand, five_prime_ext)
                    criteria.belongs_to_next_gene(peak, next_gene)
                    criteria.truncate_5_prime_end(peak, next_gene, utr)
            except criteria.CriteriaFailure as e:
                logging.debug("%s - %s" % (type(e).__name__, e))
            else:
                intersect = utr.range.intersection(map(int, sorted(spat_pileups[peak.chr], key=int))) if peak.chr in spat_pileups else None                
                colour="3"
                if intersect:                    
                    if peak.strand == "+":                        
                        utr.end = max(intersect)
                    else:
                        utr.start = min(intersect)
                    colour="4"
                if utr.is_valid():
                    logging.debug("Peak {} corresponds to 3' UTR {} of gene {}".upper().format(peak.name, utr, gene.id))
                    utr.generate_feature(gene, db, colour)
                    if peak.strand == "+":
                        gene.end = utr.end
                    else:
                        gene.start = utr.start
                    features = {"gene": gene, "utr": utr}
                    features.update({"feature_{}".format(idx): f for idx, f in enumerate(list(db.children(gene)))})
                    queue.put({gene.id: features})
                    utr_found = True
    else:
        logging.debug("No features found near peak %s" % peak.name)
        queue.put(NoNearbyFeatures)
        return
    if not utr_found:
        queue.put(None)
