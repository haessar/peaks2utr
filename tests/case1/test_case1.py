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


class TestCase1(unittest.TestCase):
    def setUp(self):
        db_path = os.path.join(TEST_DIR, "case1.db")
        gffutils.create_db(os.path.join(TEST_DIR, "case1.gtf"), db_path, force=True)
        self.db = FeatureDB(db_path)
        self.coverage_gaps = ZeroCoverageIntervalsDict()
        self.truncation_points = SPATTruncationPointsDict()
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["", ""])
        self.args.gtf_in = True

    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "case1.db"))

    def test_intronic_gene_truncation(self):
        expected_annotations = {
            "ENSMUSG00000087590": UTR(33795987, 33796107),
            "ENSMUSG00000092819": UTR(33795296, 33795382)
        }
        peaks_filename = os.path.join(TEST_DIR, "forward_peaks.broadPeak")
        peaks = BroadPeaksList(broadpeak_fn=peaks_filename, strand="forward")
        annotations = AnnotationsDict()
        pipeline = AnnotationsPipeline(peaks, self.args, queue=Queue())
        for peak in peaks:
            pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
            result = pipeline.queue.get()
            annotations.update(result)
            for gene, features in annotations.items():
                self.assertIn(gene, expected_annotations)
                self.assertEqual(features["utr"].range, expected_annotations[gene].range)


if __name__ == '__main__':
    unittest.main()
