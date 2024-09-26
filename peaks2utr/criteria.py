import logging

from .constants import FeatureTypes
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
def assert_whether_utr_already_annotated(peak, transcript, db, override_utr, extend_utr):
    """
    If the canonical annotation for this transcript already contains a 'three_prime_UTR' annotation, we either want to
    leave this as is, extend it or override it.
    """
    existing_utrs = list(db.children(transcript, featuretype=FeatureTypes.ThreePrimeUTR))
    if existing_utrs:
        if len(existing_utrs) > 1:
            logging.debug("Multiple existing 3' UTRs found for transcript %s" % transcript.id)
        if any((override_utr, extend_utr)):
            min_start = min(utr.start for utr in existing_utrs)
            max_end = max(utr.end for utr in existing_utrs)
            if transcript.strand == "+":
                transcript.end = min_start - 1 if override_utr else max_end
            else:
                transcript.start = max_end + 1 if override_utr else min_start
        else:
            raise CriteriaFailure("3' UTR already annotated for transcript %s near peak %s" % (transcript.id, peak.name))


@track_failed_peaks
def assert_peak_not_a_subset_of_transcript(peak, transcript):
    """
    If a peak occurs entirely within an existing transcript annotation (i.e. it's a subset), we consider that it is already
    accounted for and can't possibly refer to a new UTR.
    """
    if peak.range.issubset(transcript.range):
        raise CriteriaFailure("%s %s wholly contained within transcript %s"
                              % (peak.__class__.__name__, peak.name, transcript.id))


def assert_transcript_not_a_subset_of_adjacent_gene(transcript, adj_transcript, gene):
    """
    If a transcript occurs entirely within another gene's transcript, its 3' UTR should not be annotated.
    """
    if transcript.range.issubset(adj_transcript.range):
        raise CriteriaFailure("%s %s wholly contained within transcript %s of gene %s"
                              % (transcript.__class__.__name__, transcript.id, adj_transcript.id, gene.id))


@track_failed_peaks
def assert_3_prime_end_and_truncate(peak, transcript, utr):
    """
    If a peak occurs at the untranslated 3'-end of a transcript, we need to set the utr start/end to occur at the end/start of
    the existing transcript annotation, respective of strand.
    Otherwise, we take advantage of the fact that the 'assert_peak_not_a_subset_of_transcript' criteria has passed to assume
    it must correspond to the 5'-end of the transcript.
    """
    if peak.strand == "+" and peak.end > transcript.end:
        utr.start = transcript.end + 1
    elif peak.strand == "-" and peak.start < transcript.start:
        utr.end = transcript.start - 1
    else:
        raise CriteriaFailure("Peak %s corresponds to 5'-end of transcript %s" % (peak.name, transcript.id))


def truncate_to_adjacent_transcript(peak, transcript, utr, adj_transcript, gene, five_prime_ext=0):
    """
    If a peak is broad enough that it overlaps a transcript of another gene, we check for an
    intersection and truncate if it exists (taking into account assumed 5' extension).
    """
    if utr.range.intersection(adj_transcript.range):
        logging.debug("Peak %s overlapping transcript %s of gene %s: Truncating" % (peak.name, adj_transcript.id, gene.id))
        if peak.strand == "+" and adj_transcript.start > transcript.end:
            utr.end = adj_transcript.start - 1 - five_prime_ext
        elif peak.strand == "-" and adj_transcript.end < transcript.start:
            utr.start = adj_transcript.end + 1 + five_prime_ext
