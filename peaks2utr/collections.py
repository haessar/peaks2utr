import collections
from collections.abc import Sequence
import copy
import csv
import json

import gffutils

from . import constants
from .models import Peak


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
        for gid, features in self.data.items():
            for _, f in features.items():
                # gene features are redundant in GTF output
                if self.gtf_out and f.featuretype in constants.FeatureTypes.Gene:
                    continue
                formatted_f = self._apply_feature_dialect(f, gid)
                yield str(formatted_f) + '\n'
                # exon matching three_prime_UTR for GTF output
                if self.gtf_out and f.source == __package__ and f.featuretype in constants.FeatureTypes.ThreePrimeUTR:
                    exon = copy.copy(formatted_f)
                    exon.featuretype = constants.FeatureTypes.Exon[0]
                    yield str(exon) + '\n'

    @staticmethod
    def _apply_gff_dialect(feature, attrs):
        feature.dialect = constants.GFFUTILS_GFF_DIALECT
        if feature.featuretype not in constants.FeatureTypes.Gene:
            if feature.featuretype in constants.FeatureTypes.GtfTranscript:
                attrs['Parent'] = attrs.pop('gene_id')
                attrs['ID'] = attrs.pop('transcript_id')
                feature.featuretype = constants.FeatureTypes.GffTranscript[0]
            else:
                attrs.pop('gene_id')
                attrs['Parent'] = attrs.pop('transcript_id')
                attrs['ID'] = [feature.id]
        else:
            attrs['ID'] = attrs.pop('gene_id')

    @staticmethod
    def _apply_gtf_dialect(feature, attrs, gene_id=None):
        feature.dialect = constants.GFFUTILS_GTF_DIALECT
        if feature.featuretype not in constants.FeatureTypes.Gene:
            if feature.featuretype in constants.FeatureTypes.GffTranscript:
                attrs['gene_id'] = attrs.pop('Parent')
                attrs['transcript_id'] = attrs.pop('ID')
                feature.featuretype = constants.FeatureTypes.GtfTranscript[0]
            else:
                attrs['gene_id'] = gene_id
                attrs['transcript_id'] = attrs.pop('Parent')
                attrs.pop('ID')
        else:
            attrs['gene_id'] = attrs.pop('ID')

    def _apply_feature_dialect(self, feature, gene_id):
        if self.gtf_in != self.gtf_out:
            attrs = dict(feature.attributes)
            # GTF in, GFF3 out
            if not self.gtf_out and not (attrs.get('ID') and attrs.get('Parent')):
                self._apply_gff_dialect(feature, attrs)
            # GFF3 in, GTF out
            elif self.gtf_out and not (attrs.get('gene_id') and attrs.get('transcript_id')):
                self._apply_gtf_dialect(feature, attrs, gene_id)
            feature.attributes = gffutils.attributes.Attributes(**attrs)
        return feature

    def filter(self, **kwargs):
        """
        Filter features for given attributes.
        """
        all_features = [vv for v in self.values() for vv in v.values()]
        for attr, obj in kwargs.items():
            if isinstance(obj, Sequence) and not isinstance(obj, str):
                filtered_features = [f for f in all_features if getattr(f, attr) in obj]
            else:
                filtered_features = [f for f in all_features if getattr(f, attr) == obj]
        return filtered_features


class ZeroCoverageIntervalsDict(collections.UserDict):
    """
    Dictionary of zero coverage intervals per chromosome from parsed BED file.
    """
    class Interval:
        """
        1-based included
        init from 0-based half-opened
        """
        def __init__(self, start, end):
            self.start = int(start) + 1
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
            return [i for i in self[chr] if i.start <= base <= i.end]
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
    List of MACS broad peaks.
    """
    def __init__(self, initlist=None, broadpeak_fn=None, strand=None):
        super().__init__(initlist)
        if broadpeak_fn:
            with open(broadpeak_fn, 'r') as f:
                self.data = [Peak(*peak) for peak in csv.reader(f, delimiter="\t")]
                for peak in self.data:
                    peak.strand = constants.STRAND_MAP.get(strand)
