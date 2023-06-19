import os
import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import AnnotationsPipeline, NoNearbyFeatures, PotentialUTRZeroCoverage
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
        pipeline = AnnotationsPipeline(peaks, self.args, queue=Queue())
        for peak in peaks:
            if peak.name in expected_annotations:
                pipeline.annotate_utr_for_peak(self.db, peak, self.truncation_points, self.coverage_gaps)
                if expected_annotations[peak.name] is None:
                    self.assertIsNone(pipeline.queue.get())
                elif expected_annotations[peak.name] is NoNearbyFeatures:
                    self.assertIsInstance(pipeline.queue.get(), NoNearbyFeatures)
                elif expected_annotations[peak.name] is PotentialUTRZeroCoverage:
                    self.assertIsInstance(pipeline.queue.get(), PotentialUTRZeroCoverage)
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


if __name__ == '__main__':
    unittest.main()
