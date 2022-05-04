import csv
from queue import Queue
import unittest
from unittest import mock

import gffutils

from peaks2utr.annotations import Annotations, NoNearbyFeatures, annotate_utr_for_peak
from peaks2utr.models import Peak, UTR


class TestUTRAnnotation(unittest.TestCase):
    def setUp(self):
        self.db = gffutils.create_db("test/Chr1.gff", "test/Chr1.db", force=True)

    def strand_annotations(self, peaks_filename, strand, expected_annotations, max_distance):
        with open(peaks_filename, 'r') as fin:
            peaks = csv.reader(fin, delimiter="\t")
            for peak in peaks:
                peak = Peak(*peak)
                if peak.name in expected_annotations:
                    annotations = Annotations()
                    queue = Queue()
                    peak.strand = strand
                    with mock.patch('peaks2utr.annotations.cached') as cached_mock:
                        cached_mock.return_value = "test/unmapped.json"
                        annotate_utr_for_peak(self.db, queue, peak, max_distance=max_distance)
                    if expected_annotations[peak.name] is None:
                        self.assertIsNone(queue.get())                        
                    elif expected_annotations[peak.name] is NoNearbyFeatures:
                        self.assertTrue(queue.get() is NoNearbyFeatures)
                    else:
                        while not queue.empty():
                            result = queue.get()
                            if type(result) == dict:
                                annotations.update(result)
                        for gene in expected_annotations[peak.name].keys():
                            self.assertIn(gene, annotations)
                            self.assertEqual(annotations.data[gene]['utr'], expected_annotations[peak.name][gene])

    def test_forward_strand_annotations(self):
        expected_annotations = {
            'forward_peak_3': None,
            'forward_peak_6': {'PBANKA_0100041.1.1': UTR(14118, 17222)},
            'forward_peak_9': {'PBANKA_0100100.1.1': UTR(27916, 29285), 'PBANKA_0100200.1.1': UTR(30482, 31051)},
            'forward_peak_10': {'PBANKA_0100200.1.1': UTR(30482, 33095)},
        }
        peaks_filename = "test/test_forward_peaks.broadPeak"
        self.strand_annotations(peaks_filename, '+', expected_annotations, max_distance=2500)

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
        peaks_filename = "test/test_reverse_peaks.broadPeak"
        self.strand_annotations(peaks_filename, '-', expected_annotations, max_distance=2500)


if __name__ == '__main__':
    unittest.main()
