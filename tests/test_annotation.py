import os
import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import AnnotationsPipeline, NoNearbyFeatures
from peaks2utr.models import UTR, FeatureDB
from peaks2utr.collections import AnnotationsDict, BroadPeaksList, ZeroCoverageIntervalsDict, SPATTruncationPointsDict

TEST_DIR = os.path.dirname(__file__)


class TestUTRAnnotation(unittest.TestCase):
    def setUp(self):
        db_path = os.path.join(TEST_DIR, "Chr1.db")
        gffutils.create_db(os.path.join(TEST_DIR, "Chr1.gtf"), db_path, force=True)
        self.db = FeatureDB(db_path)
        self.coverage_gaps = ZeroCoverageIntervalsDict()
        self.truncation_points = SPATTruncationPointsDict()
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["", ""])
        self.args.gtf_in = True
        self.args.max_distance = 2500

    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "Chr1.db"))

    def strand_annotations(self, peaks_filename, strand, expected_annotations):
        peaks = BroadPeaksList(broadpeak_fn=peaks_filename, strand=strand)
        pipeline = AnnotationsPipeline(peaks, self.args, queue=Queue())
        for peak in peaks:
            if peak.name in expected_annotations:
                pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                if expected_annotations[peak.name] is None:
                    self.assertIsNone(pipeline.queue.get())
                elif expected_annotations[peak.name] is NoNearbyFeatures:
                    self.assertIsInstance(pipeline.queue.get(), NoNearbyFeatures)
                else:
                    result = None
                    annotations = AnnotationsDict()
                    while not pipeline.queue.empty():
                        result = pipeline.queue.get()
                        if type(result) == dict:
                            annotations.update(result)
                    for gene in expected_annotations[peak.name].keys():
                        self.assertIn(gene, annotations)
                        self.assertEqual(annotations.data[gene]['utr'].range, expected_annotations[peak.name][gene].range)

    def test_forward_strand_annotations(self):
        expected_annotations = {
            'forward_peak_3': None,
            # Up to end of forward_peak_6
            'forward_peak_6': {'PBANKA_0100041.1': UTR(14119, 17222)},
            # PBANKA_0100100.1: Up to beginning of PBANKA_0100200.1
            'forward_peak_9': {'PBANKA_0100100.1': UTR(27917, 29284), 'PBANKA_0100200.1': UTR(30483, 31051)},
            # Up to end of forward_peak_10
            'forward_peak_10': {'PBANKA_0100200.1': UTR(30483, 33095)},
        }
        peaks_filename = os.path.join(TEST_DIR, "test_forward_peaks.broadPeak")
        self.strand_annotations(peaks_filename, 'forward', expected_annotations)

    def test_reverse_strand_annotations(self):
        expected_annotations = {
            # Up to end of reverse_peak_1
            'reverse_peak_1': {'PBANKA_0100021.1': UTR(802, 1097)},
            'reverse_peak_3': NoNearbyFeatures,
            'reverse_peak_4': NoNearbyFeatures,
            # Up to end of reverse_peak_5
            'reverse_peak_5': {'PBANKA_0100061.1': UTR(21297, 21969)},
            'reverse_peak_6': None,
            # Up to end of reverse_peak_11
            'reverse_peak_11': {'PBANKA_0100800.1': UTR(45613, 47558)},
            'reverse_peak_13': None,
            'reverse_peak_14': None,
            # Up to end of reverse_peak_24
            'reverse_peak_24': {'PBANKA_0101500.1': UTR(77311, 77683)},
            # Up to gene_id PBANKA_0102100.1
            'reverse_peak_30': {'PBANKA_0102200.1': UTR(95596, 96817)},
            # Up to end of reverse_peak_54
            'reverse_peak_54': {'PBANKA_0103400.1': UTR(154887, 155565)},
            # Up to end of reverse_peak_143
            'reverse_peak_143': {'PBANKA_0111300.1': UTR(437497, 438264)}
        }
        peaks_filename = os.path.join(TEST_DIR, "test_reverse_peaks.broadPeak")
        self.strand_annotations(peaks_filename, 'reverse', expected_annotations)


if __name__ == '__main__':
    unittest.main()
