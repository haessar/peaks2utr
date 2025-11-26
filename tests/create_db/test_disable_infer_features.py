import os
import os.path
from queue import Queue
import unittest

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.annotations import AnnotationsPipeline, PotentialUTRZeroCoverage
from peaks2utr.collections import BroadPeaksList, ZeroCoverageIntervalsDict, SPATTruncationPointsDict
from peaks2utr.models import FeatureDB

TEST_DIR = os.path.dirname(__file__)


class TestDisableInferFeatures(unittest.TestCase):    
    def setUp(self):
        self.db1_path = os.path.join(TEST_DIR, "case3.no_parent.1.db")
        self.db2_path = os.path.join(TEST_DIR, "case3.no_parent.2.db")
        self.db3_path = os.path.join(TEST_DIR, "case3.no_parent.3.db")
        gffutils.create_db(os.path.join(TEST_DIR, "case3.no_parent.gtf"), self.db1_path, force=True, merge_strategy="create_unique",
                           disable_infer_genes=True,
                           disable_infer_transcripts=True)
        gffutils.create_db(os.path.join(TEST_DIR, "case3.no_parent.gtf"), self.db2_path, force=True, merge_strategy="create_unique",
                           disable_infer_genes=True,
                           disable_infer_transcripts=False)
        gffutils.create_db(os.path.join(TEST_DIR, "case3.no_parent.gtf"), self.db3_path, force=True, merge_strategy="create_unique",
                           disable_infer_genes=False,
                           disable_infer_transcripts=False)
        self.db1 = FeatureDB(self.db1_path)
        self.db2 = FeatureDB(self.db2_path)
        self.db3 = FeatureDB(self.db3_path)
        self.coverage_gaps = ZeroCoverageIntervalsDict()
        self.truncation_points = SPATTruncationPointsDict()
        forward_peaks_filename = os.path.join(TEST_DIR, "forward_peaks.broadPeak")
        self.forward_peaks = BroadPeaksList(broadpeak_fn=forward_peaks_filename, strand="forward")
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["", ""])
        self.args.gtf_in = True
    
    
    def tearDown(self):
        os.remove(self.db1_path)
        os.remove(self.db2_path)
        os.remove(self.db3_path)
    
    def _annotate_db(self, db):
        self.args.max_distance = 5000
        self.args.override_utr = True
        pipeline = AnnotationsPipeline(self.forward_peaks, self.args, queue=Queue())
        for peak in self.forward_peaks:
            pipeline.annotate_utr_for_peak(db, peak, self.truncation_points, self.coverage_gaps)
            return pipeline.queue.get()
    
    def test_utr_annotation_without_inference(self):
        # No UTR annotation when either gene or transcript is missing
        self.assertFalse(self._annotate_db(self.db1))
        self.assertFalse(self._annotate_db(self.db2))
        
        # UTR is correctly annotated when disable_infer_genes and disable_infer_transcripts are both False
        self.assertIn("ENSMUSG00000112145", self._annotate_db(self.db3))


if __name__ == '__main__':
    unittest.main()
