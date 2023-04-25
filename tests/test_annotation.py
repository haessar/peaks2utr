import os
import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import Annotations, NoNearbyFeatures
from peaks2utr.models import UTR, FeatureDB
from peaks2utr.collections import BroadPeaksList, ZeroCoverageIntervalsDict, SPATTruncationPointsDict

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
        self.args.gtf_in = False
        self.args.max_distance = 2500

    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "Chr1.db"))

    def strand_annotations(self, peaks_filename, strand, expected_annotations):
        peaks = BroadPeaksList(broadpeak_fn=peaks_filename, strand=strand)
        annotations = Annotations(peaks, self.args, queue=Queue())
        for peak in peaks:
            if peak.name in expected_annotations:
                annotations.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                if expected_annotations[peak.name] is None:
                    self.assertIsNone(annotations.queue.get())
                elif expected_annotations[peak.name] is NoNearbyFeatures:
                    self.assertTrue(annotations.queue.get() is NoNearbyFeatures)
                else:
                    result = None
                    while not annotations.queue.empty():
                        result = annotations.queue.get()
                        if type(result) == dict:
                            annotations.update(result)
                    for gene in expected_annotations[peak.name].keys():
                        self.assertIn(gene, annotations)
                        self.assertEqual(annotations.data[gene]['utr'].range, expected_annotations[peak.name][gene].range)

    def test_forward_strand_annotations(self):
        expected_annotations = {
            'forward_peak_3': None,
            'forward_peak_6': {'PBANKA_0100041.1.1': UTR(14118, 17222)},
            'forward_peak_9': {'PBANKA_0100100.1.1': UTR(27916, 29285), 'PBANKA_0100200.1.1': UTR(30482, 31051)},
            'forward_peak_10': {'PBANKA_0100200.1.1': UTR(30482, 33095)},
        }
        peaks_filename = os.path.join(TEST_DIR, "test_forward_peaks.broadPeak")
        self.strand_annotations(peaks_filename, 'forward', expected_annotations)

    def test_reverse_strand_annotations(self):
        expected_annotations = {
            'reverse_peak_1': {'PBANKA_0100021.1.1': UTR(801, 1098)},
            'reverse_peak_3': NoNearbyFeatures,
            'reverse_peak_4': NoNearbyFeatures,
            'reverse_peak_5': {'PBANKA_0100061.1.1': UTR(21296, 21970)},
            'reverse_peak_6': None,
            'reverse_peak_11': {'PBANKA_0100800.1.1': UTR(45612, 47559)},
            'reverse_peak_13': None,
            'reverse_peak_14': None,
            'reverse_peak_24': {'PBANKA_0101500.1.1': UTR(77310, 77684)},
            'reverse_peak_30': {'PBANKA_0102200.1.1': UTR(95595, 96818)},
            'reverse_peak_54': {'PBANKA_0103400.1.1': UTR(154886, 155566)},
            'reverse_peak_143': {'PBANKA_0111300.1.1': UTR(437496, 438265)}
        }
        peaks_filename = os.path.join(TEST_DIR, "test_reverse_peaks.broadPeak")
        self.strand_annotations(peaks_filename, 'reverse', expected_annotations)


if __name__ == '__main__':
    unittest.main()
