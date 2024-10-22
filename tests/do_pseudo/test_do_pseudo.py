import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import AnnotationsPipeline
from peaks2utr.collections import BroadPeaksList, ZeroCoverageIntervalsDict, SPATTruncationPointsDict
from peaks2utr.models import UTR, FeatureDB


TEST_DIR = os.path.dirname(__file__)


class TestDoPseudo(unittest.TestCase):
    def setUp(self):
        db_path = os.path.join(TEST_DIR, "PVL_12_v1.db")
        gffutils.create_db(os.path.join(TEST_DIR, "PVL_12_v1.gff"), db_path, force=True)
        self.db = FeatureDB(db_path)
        self.coverage_gaps = ZeroCoverageIntervalsDict()
        self.truncation_points = SPATTruncationPointsDict()
        forward_peaks_filename = os.path.join(TEST_DIR, "forward_peaks.broadPeak")
        self.forward_peaks = BroadPeaksList(broadpeak_fn=forward_peaks_filename, strand="forward")
        reverse_peaks_filename = os.path.join(TEST_DIR, "reverse_peaks.broadPeak")
        self.reverse_peaks = BroadPeaksList(broadpeak_fn=reverse_peaks_filename, strand="reverse")
        self.argparser = prepare_argparser()
    
    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "PVL_12_v1.db"))
    
    def test_do_pseudo_true(self):
        self.args = self.argparser.parse_args(["", "", "--do-pseudo", "--max-distance", "2000"])
        expected_annotations = {
            "forward_peak_3": {"PVL_120054200": UTR(2520543, 2522849)},
            "reverse_peak_3": {"PVL_120054300": UTR(2522147, 2522539)}
        }
        self._run_test(expected_annotations)
    
    def test_do_pseudo_false(self):
        self.args = self.argparser.parse_args(["", "", "--max-distance", "2000"])
        expected_annotations = {
            "forward_peak_3": None,
            "reverse_peak_3": None
        }
        self._run_test(expected_annotations)
    
    def test_do_pseudo_true_no_strand_overlap_annotations(self):
        self.args = self.argparser.parse_args(["", "", "--do-pseudo", "--max-distance", "2000", "--no-strand-overlap"])
        expected_annotations = {
            "forward_peak_3": {"PVL_120054200": UTR(2520543, 2522539)},
            "reverse_peak_3": {"PVL_120054300": UTR(2522147, 2522539)}
        }
        self._run_test(expected_annotations)

    def _run_test(self, expected_annotations):
        pipeline = AnnotationsPipeline(self.forward_peaks + self.reverse_peaks, self.args, queue=Queue())
        for peak in self.forward_peaks + self.reverse_peaks:
            if peak.name in expected_annotations:
                pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                result = pipeline.queue.get()
                if result:
                    for gene in expected_annotations[peak.name].keys():
                        self.assertIn(gene, result)
                        self.assertEqual(expected_annotations[peak.name][gene].range, result[gene]["utr"].range)
                else:
                    assert expected_annotations[peak.name] is None
        

if __name__ == '__main__':
    unittest.main()
