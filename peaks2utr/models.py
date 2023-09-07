"""
For a discussion of 0-based/1-based counting systems,
see https://genome-blog.soe.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
"""
import re

import gffutils

from .constants import AnnotationColour, FeatureTypes, STRAND_CIGAR_SOFT_CLIP_REGEX, GFFUTILS_GFF_DIALECT, GFFUTILS_GTF_DIALECT


class RangeMixin:
    """
    Like gff/gtf this class mixin is 1-based included
    """
    start: int
    end: int

    @property
    def range(self):
        return set(range(self.start, self.end + 1))

    @property
    def length(self):
        return self.end - self.start + 1


class Peak(RangeMixin):
    """
    MACS peak in BED6+3 format but 1-based included
    init from 0-based half-opened BED6+3 arguments.
    """
    def __init__(self, *args):
        self.chr = str(args[0])
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        self.name = str(args[3])
        self.score = int(args[4])
        self.strand = str(args[5])
        self.signalValue = float(args[6])
        self.pValue = float(args[7])
        self.qValue = float(args[8])

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__, str(self.__dict__))


class Feature(gffutils.Feature, RangeMixin):
    """
    gffutils.Feature with range property
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.keep_order = True


class FeatureDB(gffutils.FeatureDB):
    def _feature_returner(self, **kwargs):
        """
        Overwrite the gffutils.FeatureDB._feature_returner method to "slot in" Feature with added range property
        """
        kwargs.setdefault('dialect', self.dialect)
        kwargs.setdefault('keep_order', self.keep_order)
        kwargs.setdefault('sort_attribute_values', self.sort_attribute_values)
        return Feature(**kwargs)


class UTR(RangeMixin):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.feature = None

    def __str__(self):
        return str(self.feature) if self.feature else super().__str__()

    def __repr__(self):
        if self.feature:
            return self.feature.__repr__()
        return "<%s: (%s, %s)>" % (self.__class__.__name__, self.start, self.end)

    def __eq__(self, other):
        return self.range == other.range

    def _create_id(self, transcript, db):
        existing_utrs = list(db.children(transcript, featuretype=FeatureTypes.ThreePrimeUTR)) + \
                        list(db.children(transcript, featuretype=FeatureTypes.FivePrimeUTR))
        if existing_utrs:
            max_utr = sorted([utr.id for utr in existing_utrs], reverse=True)[0]
            max_idx = int(max_utr[-1])
            max_utr_basename = max_utr[:-1]
            return max_utr_basename + str(max_idx + 1)
        else:
            return "utr_" + transcript.id + "_1"

    def generate_feature(self, gene, transcript, db, colour=AnnotationColour.Extended, gtf_in=False):
        """
        Generate three_prime_UTR feature in gff3 format.
        """
        d = {
            "seqid": transcript.chrom,
            "source": __package__,
            "featuretype": FeatureTypes.ThreePrimeUTR[0],
            "start": self.start,
            "end": self.end,
            "score": '.',
            "strand": transcript.strand,
            "frame": '.',
            "dialect": GFFUTILS_GTF_DIALECT if gtf_in else GFFUTILS_GFF_DIALECT,
        }
        attrs = {}
        id = self._create_id(transcript, db)
        if gtf_in:
            attrs["gene_id"] = [gene.id]
            attrs["transcript_id"] = [transcript.id]
        else:
            attrs["ID"] = [id]
            attrs["Parent"] = [transcript.id]
        attrs.update({'colour': [colour]})
        d.update({"attributes": attrs})

        self.feature = Feature(id=id, **d)

    def is_valid(self):
        return self.end >= self.start


class SoftClippedRead:
    """
    Read in SAM file format and store in 1-based included
    init from 0-based half-opened
    """
    def __init__(self, chr, start, end, cigar, seq, strand):
        self.chr = str(chr)
        self.start = int(start) + 1
        self.end = int(end)
        self.cigar = str(cigar)
        self.seq = str(seq)
        self.strand = strand

    @property
    def len_soft_clipped(self):
        """
        Length of soft-clipped bases at end of read.
        """
        pattern = STRAND_CIGAR_SOFT_CLIP_REGEX.get(self.strand)
        matches = re.search(pattern, self.cigar)
        if matches:
            return int(matches.group(1))
        else:
            return 0

    @property
    def extremity(self):
        """
        Furthest base from transcript, accounting for strand.
        """
        if self.strand == "reverse":
            return self.start
        return self.end

    def poly_tail_exists(self, tail_len=10):
        """
        Return True if a poly-A/T tail of tail_len bases exists in soft-clipped portion of read.
        """
        if self.len_soft_clipped > 0:
            soft_clipped = self.seq[:self.len_soft_clipped] if self.strand == "reverse" else self.seq[-self.len_soft_clipped:]
            if "T"*tail_len in soft_clipped or "A"*tail_len in soft_clipped:
                return True
        return False
