import logging

from .utils import Counter

class CriteriaFailure(Exception):
    pass


def track_failed_peaks(f):
    """
    Decorator to track set of peaks that fail this criterion.
    """
    def wrapped(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except CriteriaFailure:
            peak = kwargs.get('peak', args[0])
            wrapped.fails.add(peak.name)
            raise
    wrapped.fails = Counter()
    return wrapped


@track_failed_peaks
def assert_whether_utr_already_annotated(peak, gene, db, override_utr, extend_utr):
    """
    If the canonical annotation for this gene already contains a 'three_prime_UTR' annotation, we either want to
    leave this as is, extend it or override it.
    """
    existing_utrs = list(db.children(gene, featuretype=['three_prime_UTR', 'three_prime_utr']))
    if existing_utrs:
        if len(existing_utrs) > 1:
            logging.debug("Multiple existing 3' UTRs found for gene %s" % gene.id) 
        if any((override_utr, extend_utr)):
            min_start = min(utr.start for utr in existing_utrs)
            min_end = min(utr.end for utr in existing_utrs)
            max_start = max(utr.start for utr in existing_utrs)
            max_end = max(utr.end for utr in existing_utrs)
            if gene.strand == "+":
                gene.end = max_start if override_utr else max_end            
            else:
                gene.start = min_end if override_utr else min_start            
        else:
            raise CriteriaFailure("3' UTR already annotated for gene %s near peak %s" % (gene.id, peak.name))


@track_failed_peaks
def assert_not_a_subset(peak, gene):
    """
    If a peak occurs entirely within an existing gene annotation (i.e. it's a subset), we consider that it is already
    accounted for and can't possibly refer to a new UTR.
    """
    if peak.range.issubset(gene.range):
        raise CriteriaFailure("Peak %s wholly contained within gene %s" % (peak.name, gene.id))


@track_failed_peaks
def assert_3_prime_end_and_truncate(peak, gene, utr):
    """
    If a peak occurs at the untranslated 3'-end of a gene, we need to set the utr start/end to occur at the end/start of
    the existing gene annotation, respective of strand.
    Otherwise, we take advantage of the fact that the 'assert_not_a_subset' criteria has passed to assume it must
    correspond to the 5'-end of the gene.
    """
    if peak.strand == "+" and peak.end > gene.end:
        utr.start = gene.end
    elif peak.strand == "-" and peak.start < gene.start:
        utr.end = gene.start
    else:
        raise CriteriaFailure("Peak %s corresponds to 5'-end of gene %s" % (peak.name, gene.id))


def truncate_5_prime_end(peak, next_gene, utr):
    """
    If a peak is broad enough it can potentially overlap the 5'-end of the following gene, so we check for an
    intersection and truncate if it exists.
    """
    if utr.range.intersection(next_gene.range):
        logging.debug("Peak %s overlapping following gene %s: Truncating" % (peak.name, next_gene.id))
        if peak.strand == "+":
            utr.end = next_gene.start
        else:
            utr.start = next_gene.end


def belongs_to_next_gene(peak, next_gene):
    """
    If the max_distance is large enough, it's entirely possible for a peak to occur after the start of the following
    gene and still be within range of the present gene. In this case we want to consider it "belonging" to the following
    gene only.
    """
    if (peak.strand == "+" and peak.start > next_gene.start) or (peak.strand == "-" and peak.end < next_gene.end):
        raise CriteriaFailure("Peak %s belongs entirely to following gene %s" % (peak.name, next_gene.id))
