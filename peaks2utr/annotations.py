import collections
import copy
import logging

from . import criteria, models


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


def set_gene_range(gene, strand, five_prime_ext=0):
    """Assign a range to gene taking into account assumed 5' extension"""
    if strand == "+":
        gene.start -= five_prime_ext
    else:
        gene.end += five_prime_ext
    gene.range = range(gene.start, gene.end)


def annotate_utr_for_peak(db, annotations, peak, max_distance, override_utr, five_prime_ext):
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
                gene = copy.deepcopy(genes[idx])
                set_gene_range(gene, peak.strand)
                if not override_utr:
                    criteria.assert_not_already_annotated(peak, gene, db)
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
                if utr.is_valid():
                    logging.debug("Peak {} corresponds to 3' UTR {} of gene {}".upper().format(peak.name, utr, gene.id))
                    utr.generate_feature(gene)
                    annotations[gene.id] = utr
    else:
        logging.debug("No features found near peak %s" % peak.name)
        annotations.no_features_counter += 1
