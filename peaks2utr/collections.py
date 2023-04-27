import collections
import csv
import json

import gffutils

from . import constants
from .models import Peak
from .utils import feature_from_line


class AnnotationsDict(collections.UserDict):
    """
    Dictionary of features per gene id.
    """
    def __init__(self, dict=None, args=None):
        super().__init__(dict)
        if args:
            self.gtf_in = args.gtf_in
            self.gtf_out = args.gtf_out

    def __setitem__(self, gene, new_features):
        existing_features = self.get(gene)
        if existing_features:
            if new_features["utr"].range.issubset(existing_features["utr"].range):
                return
        self.data[gene] = new_features

    def iter_feature_strings(self):
        for _, features in self.data.items():
            yield '\n'.join([str(self._apply_feature_dialect(f)) for _, f in features.items()]) + '\n'

    def _apply_feature_dialect(self, feature):
        if self.gtf_in != self.gtf_out:
            attrs = dict(feature.attributes)
            # GTF in, GFF3 out
            if not self.gtf_out and not (attrs.get('ID') and attrs.get('Parent')):
                if not attrs.get('ID'):
                    attrs['ID'] = [feature.id]
                if not attrs.get('Parent') and feature.featuretype not in constants.FeatureTypes.Gene:
                    attrs['Parent'] = attrs['gene_id'] if feature.featuretype in constants.FeatureTypes.Transcript else \
                                      attrs['transcript_id']
            # GFF3 in, GTF out
            elif self.gtf_out and not (attrs.get('gene_id') and attrs.get('transcript_id')):
                if feature.featuretype not in constants.FeatureTypes.Gene:
                    if feature.featuretype in constants.FeatureTypes.Transcript:
                        attrs['gene_id'] = attrs['Parent']
                        attrs['transcript_id'] = [feature.id]
                    else:
                        attrs['transcript_id'] = attrs['Parent']
                else:
                    attrs['gene_id'] = [feature.id]

            feature.attributes = gffutils.attributes.Attributes(**attrs)
        return feature_from_line(
            str(feature),
            dialect_in=constants.GFFUTILS_GTF_DIALECT if self.gtf_in else constants.GFFUTILS_GFF_DIALECT,
            dialect_out=constants.GFFUTILS_GTF_DIALECT if self.gtf_out else constants.GFFUTILS_GFF_DIALECT,
        )


class ZeroCoverageIntervalsDict(collections.UserDict):
    """
    Dictionary of zero coverage intervals per chromosome from parsed BED file.
    """
    class Interval:
        def __init__(self, start, end):
            self.start = int(start)
            self.end = int(end)

    def __init__(self, dict=None, bed_fn=None):
        super().__init__(dict)
        if bed_fn:
            with open(bed_fn, 'r') as f:
                for line in f.readlines():
                    chr, start, end = line.strip().split('\t')
                    if chr not in self.data:
                        self.data[chr] = []
                    self.data[chr].append(self.Interval(start, end))

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


class SPATTruncationPointsDict(collections.UserDict):
    """
    Dictionary of SPAT "truncation points" per chromosome from json file.
    """
    def __init__(self, dict=None, json_fn=None):
        super().__init__(dict)
        if json_fn:
            with open(json_fn, 'r') as f:
                self.data.update(json.load(f) or {})


class BroadPeaksList(collections.UserList):
    """
    List of MACS3 broad peaks
    """
    def __init__(self, initlist=None, broadpeak_fn=None, strand=None):
        super().__init__(initlist)
        if broadpeak_fn:
            with open(broadpeak_fn, 'r') as f:
                self.data = [Peak(*peak) for peak in csv.reader(f, delimiter="\t")]
                for peak in self.data:
                    peak.strand = constants.STRAND_MAP.get(strand)
