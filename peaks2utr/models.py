from abc import ABC
import collections
import json
import re

import gffutils

from .constants import STRAND_CIGAR_SOFT_CLIP_REGEX


class RangeMixin(ABC):
    start: int
    end: int

    @property
    def range(self):
        return set(range(self.start, self.end))


class Peak(RangeMixin):
    """
    MACS3 peak in BED 6+3 format.
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

    def _get_id(self, gene, db):
        existing_utrs = list(db.children(gene, featuretype=['three_prime_UTR', 'three_prime_utr'])) + list(db.children(gene, featuretype=['five_prime_UTR', 'five_prime_utr']))
        if existing_utrs:
            max_utr = sorted([utr.id for utr in existing_utrs], reverse=True)[0]
            max_idx = int(max_utr[-1])
            max_utr_basename = max_utr[:-1]
            return max_utr_basename + str(max_idx + 1)
        else:
            return "utr_" + gene.id + ":mRNA_1"

    def generate_feature(self, gene, db, colour="3"):
        """
        Generate three_prime_UTR feature in gff3 format.
        """
        attrs = dict(gene.attributes)
        attrs.pop('ID', None)
        attrs["ID"] = [self._get_id(gene, db)]
        attrs['Parent'] = [gene.id]
        attrs['colour'] = [colour]
        self.feature = gffutils.Feature(
            seqid=gene.chrom,
            source=__package__,
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


class SoftClippedRead:
    """
    Read in SAM file format.
    """
    def __init__(self, chr, start, end, cigar, seq, strand):
        self.chr = str(chr)
        self.start = int(start)
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

class ZeroCoverageIntervals(collections.UserDict):
    """
    Dictionary of zero coverage intervals per chromosome from parsed BED file.
    """
    class Interval:
        def __init__(self, start, end):
            self.start = int(start)
            self.end = int(end)

    def __init__(self, bed_fn=None):
        super().__init__()
        if bed_fn:
            with open(bed_fn, 'r') as f:
                for l in f.readlines():
                    chr, start, stop = l.strip().split('\t')
                    if chr not in self.data:
                        self.data[chr] = []
                    self.data[chr].append(self.Interval(start, stop))

    def filter(self, chr, base):
        """
        Filter intervals that contain base.

        This was motivated by https://github.com/haessar/peaks2utr/issues/9 - the native
        pybedtools.BedTool.filter method was causing "Too many files open" errors in some
        distributed systems.
        """
        if chr in self:
            return [i for i in self[chr] if i.start < base < i.end]
        return []

class SPATTruncationPoints(collections.UserDict):
    """
    Dictionary of SPAT "truncation points" per chromosome from json file.
    """
    def __init__(self, json_fn=None):
        super().__init__()
        if json_fn:
            with open(json_fn, 'r') as f:
                self.data.update(json.load(f))
