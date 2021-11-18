import csv
import unittest

import gffutils

from peaks_to_3pUTR.peaks_to_UTR import Annotations, Peak, UTR, annotate_utr_for_peak


class TestUTRAnnotation(unittest.TestCase):
    def setUp(self):
        self.db = gffutils.create_db("Chr1.gff", "Chr1.db", force=True)

    def strand_annotations(self, peaks_filename, strand, expected_annotations, max_distance):
        with open(peaks_filename, 'r') as fin:
            peaks = csv.reader(fin, delimiter="\t")
            for peak in peaks:
                peak = Peak(*peak)
                if peak.name in expected_annotations:
                    annotations = Annotations()
                    peak.strand = strand
                    annotate_utr_for_peak(self.db, annotations, peak, max_distance=max_distance)
                    if not expected_annotations[peak.name]:
                        self.assertFalse(bool(annotations))
                    else:
                        for gene in expected_annotations[peak.name].keys():
                            self.assertIn(gene, annotations)
                        self.assertDictEqual(annotations.data, expected_annotations[peak.name])

    def test_forward_strand_annotations(self):
        expected_annotations = {
            'forward_peak_3': {},
            'forward_peak_6': {'PBANKA_0100041.1.1': UTR(14118, 17222)},
            'forward_peak_9': {'PBANKA_0100100.1.1': UTR(27916, 29285), 'PBANKA_0100200.1.1': UTR(30482, 31051)},
            'forward_peak_10': {'PBANKA_0100200.1.1': UTR(30482, 33095)},
        }
        peaks_filename = "test_forward_peaks.broadPeak"
        self.strand_annotations(peaks_filename, '+', expected_annotations, max_distance=2500)

    def test_reverse_strand_annotations(self):
        expected_annotations = {
            'reverse_peak_1': {'PBANKA_0100021.1.1': UTR(801, 1098)},
            'reverse_peak_3': {},
            'reverse_peak_4': {},
            'reverse_peak_5': {'PBANKA_0100061.1.1': UTR(21296, 21970)},
            'reverse_peak_6': {},
            'reverse_peak_11': {'PBANKA_0100800.1.1': UTR(45612, 47559)},
            'reverse_peak_13': {},
            'reverse_peak_14': {},
            'reverse_peak_24': {'PBANKA_0101500.1.1': UTR(77310, 77684)},
            'reverse_peak_30': {'PBANKA_0102200.1.1': UTR(95595, 96818)},
            'reverse_peak_54': {'PBANKA_0103400.1.1': UTR(154886, 155566)},
            'reverse_peak_143': {'PBANKA_0111300.1.1': UTR(437496, 438265)}
        }
        peaks_filename = "test_reverse_peaks.broadPeak"
        self.strand_annotations(peaks_filename, '-', expected_annotations, max_distance=2500)


if __name__ == '__main__':
    unittest.main()
