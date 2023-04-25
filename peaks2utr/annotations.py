import collections
import copy
import logging
import math
import multiprocessing
import sqlite3

import gffutils
from tqdm import tqdm

from . import constants, criteria
from .constants import AnnotationColour, STRAND_MAP
from .collections import SPATTruncationPointsDict, ZeroCoverageIntervalsDict
from .models import UTR, FeatureDB
from .utils import cached, feature_from_line, iter_batches


class NoNearbyFeatures:
    pass


class Annotations(collections.UserDict):
    def __init__(self, peaks, args, queue=None):
        super().__init__()
        self.no_features_counter = 0
        self.peaks = peaks
        self.total_peaks = len(peaks)
        self.args = args
        self.queue = queue or multiprocessing.Queue()

    def __setitem__(self, gene, new_features):
        existing_features = self.get(gene)
        if existing_features:
            if new_features["utr"].range.issubset(existing_features["utr"].range):
                return
        self.data[gene] = new_features

    def __iter__(self):
        for _, features in self.data.items():
            yield '\n'.join([str(self._apply_feature_dialect(f)) for _, f in features.items()]) + '\n'

    def _apply_feature_dialect(self, feature):
        if self.args.gtf_in != self.args.gtf_out:
            attrs = dict(feature.attributes)
            # GTF in, GFF3 out
            if not self.args.gtf_out and not (attrs.get('ID') and attrs.get('Parent')):
                if not attrs.get('ID'):
                    attrs['ID'] = [feature.id]
                if not attrs.get('Parent') and feature.featuretype not in constants.FeatureTypes.Gene:
                    attrs['Parent'] = attrs['gene_id'] if feature.featuretype in constants.FeatureTypes.Transcript else \
                                      attrs['transcript_id']
            # GFF3 in, GTF out
            elif self.args.gtf_out and not (attrs.get('gene_id') and attrs.get('transcript_id')):
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
            dialect_in=constants.GFFUTILS_GTF_DIALECT if self.args.gtf_in else constants.GFFUTILS_GFF_DIALECT,
            dialect_out=constants.GFFUTILS_GTF_DIALECT if self.args.gtf_out else constants.GFFUTILS_GFF_DIALECT,
        )

    def __call__(self, db_path):
        """
        Args:
            db_path (string): GFFUtils sqlite3 db file path.
        """
        self.db_path = db_path

    def __enter__(self):
        self.processes = [
            self._batch_annotate_strand(batch)
            for batch in iter_batches(self.peaks, math.ceil(self.total_peaks/self.args.processors))
        ]
        for p in self.processes:
            p.start()
        self.pbar = tqdm(total=self.total_peaks, desc=f'{"INFO": <8} Iterating over peaks to annotate 3\' UTRs.')
        return self

    def __exit__(self, type, value, traceback):
        self.pbar.close()

    def _connect_db(self):
        db = sqlite3.connect(self.db_path, check_same_thread=False)
        return FeatureDB(db)

    def _batch_annotate_strand(self, peaks_batch):
        """
        Create multiprocessing Process to handle batch of peaks. Connect to sqlite3 db for each batch to prevent
        serialization issues.
        """
        truncation_points = {}
        coverage_gaps = {}
        for strand, symbol in STRAND_MAP.items():
            truncation_points[symbol] = SPATTruncationPointsDict(json_fn=cached(strand + "_unmapped.json"))
            coverage_gaps[symbol] = ZeroCoverageIntervalsDict(bed_fn=cached(strand + "_coverage_gaps.bed"))
        db = self._connect_db()
        return multiprocessing.Process(target=self._iter_peaks, args=(db, peaks_batch, truncation_points, coverage_gaps))

    def _iter_peaks(self, db, peaks_batch, truncation_points, coverage_gaps):
        for peak in peaks_batch:
            self.annotate_utr_for_peak(
                db,
                peak,
                truncation_points.get(peak.strand),
                coverage_gaps.get(peak.strand))

    def _filter_db(self, db, chr, start, end, strand, featuretype):
        features = list(db.region(
            seqid=chr,
            start=start - self.args.max_distance,
            end=end + self.args.max_distance,
            strand=strand,
            featuretype=featuretype)
        )
        return sorted(features, key=lambda x: x.start, reverse=False if strand == "+" else True)

    def annotate_utr_for_peak(self, db, peak, truncation_points, coverage_gaps):
        """
        Find genes in region of given peak and apply criteria to determine if 3' UTR exists for each.
        If so, add to multiprocessing Queue.

        Args:
            db (gffutils.interface.FeatureDB)
        """
        utr_found = False
        genes = self._filter_db(db, peak.chr, peak.start, peak.end, peak.strand, constants.FeatureTypes.Gene) or []
        if genes:
            for idx, gene in enumerate(genes):
                transcripts = db.children(
                    gene,
                    featuretype=constants.FeatureTypes.Transcript,
                    order_by="end" if peak.strand == "+" else "start",
                    reverse=True if peak.strand == "+" else False
                )
                # Take outermost transcript
                try:
                    transcript = next(transcripts)
                except StopIteration:
                    continue
                try:
                    criteria.assert_whether_utr_already_annotated(peak, transcript, db,
                                                                  self.args.override_utr, self.args.extend_utr)
                    criteria.assert_not_a_subset(peak, transcript)
                    utr = UTR(start=peak.start, end=peak.end)
                    criteria.assert_3_prime_end_and_truncate(peak, transcript, utr)
                    if len(genes) > idx + 1:
                        next_gene = copy.deepcopy(genes[idx + 1])
                        criteria.belongs_to_next_gene(peak, next_gene, self.args.five_prime_ext)
                        criteria.truncate_5_prime_end(peak, next_gene, utr, self.args.five_prime_ext)
                except criteria.CriteriaFailure as e:
                    logging.debug("%s - %s" % (type(e).__name__, e))
                else:
                    colour = AnnotationColour.Extended
                    intersect = utr.range.intersection(map(int, sorted(truncation_points[peak.chr], key=int))) \
                        if peak.chr in truncation_points else None
                    if peak.strand == "+":
                        gaps = coverage_gaps.filter(peak.chr, utr.end)
                        try:
                            gap_edge = min([g.start for g in gaps])
                        except ValueError:
                            pass
                        else:
                            utr.end = max(transcript.end, gap_edge)
                            colour = AnnotationColour.TruncatedZeroCoverage
                    else:
                        gaps = coverage_gaps.filter(peak.chr, utr.start)
                        try:
                            gap_edge = max([g.end for g in gaps])
                        except ValueError:
                            pass
                        else:
                            utr.start = min(transcript.start, gap_edge)
                            colour = AnnotationColour.TruncatedZeroCoverage
                    if intersect:
                        if peak.strand == "+":
                            utr.end = max(intersect)
                        else:
                            utr.start = min(intersect)
                        colour = AnnotationColour.ExtendedWithSPAT
                    if utr.is_valid():
                        logging.debug("Peak {} corresponds to 3' UTR {} of gene {}".upper().format(peak.name, utr, gene.id))
                        utr.generate_feature(gene, transcript, db, colour, self.args.gtf_in, self.args.gtf_out)
                        features = {"gene": gene, "transcript": transcript}
                        features.update({"feature_{}".format(idx): f for idx, f in enumerate(list(db.children(transcript)))
                                        if f.id != transcript.id and f.id != gene.id})
                        features.update({"utr": utr.feature})
                        if peak.strand == "+":
                            gene.end = transcript.end = utr.end
                        else:
                            gene.start = transcript.start = utr.start
                        self.queue.put({transcript.id: features})
                        utr_found = True
        else:
            logging.debug("No features found near peak %s" % peak.name)
            self.queue.put(NoNearbyFeatures)
            return
        if not utr_found:
            self.queue.put(None)
