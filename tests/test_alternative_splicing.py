import os
import os.path
import unittest

from asgiref.sync import async_to_sync

from peaks2utr import prepare_argparser
from peaks2utr.constants import FeatureTypes, GFFUTILS_GTF_DIALECT
from peaks2utr.models import Feature, FeatureDB
from peaks2utr.preprocess import create_db

TEST_DIR = os.path.dirname(__file__)


class TestAlternativeSplicing(unittest.TestCase):
    def setUp(self):
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["", ""])
        self.args.alternative_splicing = False
        chr = "ateAlb2.1"
        strand = "+"
        gene1_id = 'Aalb_00004275'
        gene2_id = 'Aalb_00004276'
        gene3_id = 'Aalb_00004277'
        self.transcript1 = Feature(chr, id=gene1_id + '-T1', featuretype=FeatureTypes.GtfTranscript[0],
                                   start=153660, end=157786, strand=strand, dialect=GFFUTILS_GTF_DIALECT,
                                   attributes={"transcript_id": [gene1_id + '-T1'], "gene_id": [gene1_id]})
        self.transcript2 = Feature(chr, id=gene2_id + '-T1', featuretype=FeatureTypes.GtfTranscript[0],
                                   start=153660, end=162595, strand=strand, dialect=GFFUTILS_GTF_DIALECT,
                                   attributes={"transcript_id": [gene2_id + '-T1'], "gene_id": [gene2_id]})
        self.transcript3 = Feature(chr, id=gene3_id + '-T1', featuretype=FeatureTypes.GtfTranscript[0],
                                   start=157165, end=441035, strand=strand, dialect=GFFUTILS_GTF_DIALECT,
                                   attributes={"transcript_id": [gene3_id + '-T1'], "gene_id": [gene3_id]})
        self.db_path = os.path.join(TEST_DIR, "alternative_splicing.db")
        async_to_sync(create_db)(os.path.join(TEST_DIR, "alternative_splicing.gtf"), self.db_path, alternative_splicing=False)
        self.db = FeatureDB(self.db_path)

    def tearDown(self):
        os.remove(self.db_path)


if __name__ == '__main__':
    unittest.main()
