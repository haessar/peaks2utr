import os
import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import AnnotationsPipeline
from peaks2utr.collections import BroadPeaksList, ZeroCoverageIntervalsDict, SPATTruncationPointsDict
from peaks2utr.models import UTR, FeatureDB

TEST_DIR = os.path.dirname(__file__)


class TestCase2(unittest.TestCase):
    def setUp(self):
        db_path = os.path.join(TEST_DIR, "case2.db")
        gffutils.create_db(os.path.join(TEST_DIR, "case2.gtf"), db_path, force=True)
        self.db = FeatureDB(db_path)
        self.coverage_gaps = ZeroCoverageIntervalsDict()
        self.truncation_points = SPATTruncationPointsDict()
        forward_peaks_filename = os.path.join(TEST_DIR, "forward_peaks.broadPeak")
        self.forward_peaks = BroadPeaksList(broadpeak_fn=forward_peaks_filename, strand="forward")
        reverse_peaks_filename = os.path.join(TEST_DIR, "reverse_peaks.broadPeak")
        self.reverse_peaks = BroadPeaksList(broadpeak_fn=reverse_peaks_filename, strand="reverse")
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["", ""])
        self.args.gtf_in = True

    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "case2.db"))

    def test_gene_within_exon(self):
        expected_annotations = {"forward_peak_32164": None}
        pipeline = AnnotationsPipeline(self.forward_peaks, self.args, queue=Queue())
        for peak in self.forward_peaks:
            if peak.name in expected_annotations:
                pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                result = pipeline.queue.get()
                self.assertEqual(expected_annotations[peak.name], result)

    def test_override_utr(self):
        self.args.max_distance = 5000
        self.args.override_utr = True
        expected_annotations = {
            "reverse_peak_32037": {"ENSMUSG00000033396": UTR(122048355, 122053879)},
            "forward_peak_32170": {"ENSMUSG00000027236": UTR(122052099, 122053546)}
        }
        pipeline = AnnotationsPipeline(self.forward_peaks, self.args, queue=Queue())
        for peak in self.forward_peaks + self.reverse_peaks:
            if peak.name in expected_annotations:
                pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                result = pipeline.queue.get()
                if result:
                    gene = [g for g in expected_annotations[peak.name].keys()][0]
                    self.assertIn(gene, result)
                    self.assertEqual(expected_annotations[peak.name][gene].range, result[gene]["utr"].range)


if __name__ == '__main__':
    unittest.main()
