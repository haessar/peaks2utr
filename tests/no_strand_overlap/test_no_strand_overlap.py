import os
import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import AnnotationsPipeline
from peaks2utr.collections import AnnotationsDict, BroadPeaksList, ZeroCoverageIntervalsDict, SPATTruncationPointsDict
from peaks2utr.models import UTR, FeatureDB

TEST_DIR = os.path.dirname(__file__)


class TestNoStrandOverlap(unittest.TestCase):
    def setUp(self):
        db_path = os.path.join(TEST_DIR, "Chr1.db")
        gffutils.create_db(os.path.join(TEST_DIR, "Chr1.gff"), db_path, force=True, disable_infer_genes=True, disable_infer_transcripts=True)
        self.db = FeatureDB(db_path)
        self.coverage_gaps = ZeroCoverageIntervalsDict()
        self.truncation_points = SPATTruncationPointsDict()
        forward_peaks_filename = os.path.join(TEST_DIR, "forward_peaks.broadPeak")
        self.forward_peaks = BroadPeaksList(broadpeak_fn=forward_peaks_filename, strand="forward")
        reverse_peaks_filename = os.path.join(TEST_DIR, "reverse_peaks.broadPeak")
        self.reverse_peaks = BroadPeaksList(broadpeak_fn=reverse_peaks_filename, strand="reverse")
        self.argparser = prepare_argparser()

    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "Chr1.db"))

    def _do_annotation(self, args, expected_annotations):
        pipeline = AnnotationsPipeline(self.forward_peaks + self.reverse_peaks, args, queue=Queue())
        for peak in self.forward_peaks + self.reverse_peaks:
            if peak.name in expected_annotations:
                pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                result = None
                annotations = AnnotationsDict()
                while not pipeline.queue.empty():
                    result = pipeline.queue.get()
                    if type(result) == dict:
                        annotations.update(result)
                for gene in expected_annotations[peak.name].keys():
                    self.assertIn(gene, annotations)
                    self.assertEqual(annotations.data[gene]['utr'].range, expected_annotations[peak.name][gene].range)

    def test_utr_no_overlap(self):
        args = self.argparser.parse_args(["", "", "--max-distance", "2000", "--no-strand-overlap"])
        expected_annotations = {
            "forward_peak_104": {"C4B63_1g127": UTR(278110, 280146)},
            "reverse_peak_99": {"C4B63_1g1319": UTR(279445, 280146)}
        }
        self._do_annotation(args, expected_annotations)
    
    def test_utr_no_overlap_do_pseudo(self):
        args = self.argparser.parse_args(["", "", "--max-distance", "2000", "--no-strand-overlap", "--do-pseudo"])
        expected_annotations = {
            "forward_peak_104": {"C4B63_1g127": UTR(278110, 280028)},
            "reverse_peak_99": {
                "C4B63_1g1319": UTR(279445, 280146),
                "C4B63_1g1320": UTR(279445, 280028)}
        }
        self._do_annotation(args, expected_annotations)
    
    def test_utr_with_overlap(self):
        args = self.argparser.parse_args(["", "", "--max-distance", "2000"])
        expected_annotations = {
            "forward_peak_104": {"C4B63_1g127": UTR(278110, 281990)},
            "reverse_peak_99": {"C4B63_1g1319": UTR(279445, 280146)}
        }
        self._do_annotation(args, expected_annotations)

if __name__ == '__main__':
    unittest.main()
