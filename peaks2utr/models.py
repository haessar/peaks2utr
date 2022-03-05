from abc import ABC
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
        """
        Generate three_prime_UTR feature in gff3 format.
        """
        attrs = dict(gene.attributes)
        attrs.pop('ID', None)
        attrs['Parent'] = [gene.id]
        attrs['colour'] = ['3']
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
    Read in SAM file format
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
        pattern = STRAND_CIGAR_SOFT_CLIP_REGEX.get(self.strand)
        matches = re.search(pattern, self.cigar)
        if matches:
            return int(matches.group(1))
        else:
            return 0
    
    @property
    def extremity(self):
        if self.strand == "reverse":
            return self.start
        return self.end

    def poly_tail_exists(self, tail_len=10):
        if self.len_soft_clipped > 0:
            soft_clipped = self.seq[:self.len_soft_clipped] if self.strand == "reverse" else self.seq[-self.len_soft_clipped:]
            if "T"*tail_len in soft_clipped or "A"*tail_len in soft_clipped:                
                return True
        return False
