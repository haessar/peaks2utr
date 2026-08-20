import gzip
import os
import os.path
import shutil
import tempfile
import unittest

import gffutils

from peaks2utr.models import FeatureDB

TEST_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.join(TEST_DIR, "..")


class TestCreateDbWithCompressedInput(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.mkdtemp()

        # Gzip the GTF fixture
        cls.gtf_gz = os.path.join(cls.temp_dir, "case3.no_parent.gtf.gz")
        with open(os.path.join(TEST_DIR, "case3.no_parent.gtf"), "rb") as f_in:
            with gzip.open(cls.gtf_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Gzip the GFF3 fixture
        cls.gff_gz = os.path.join(cls.temp_dir, "Chr1.gff.gz")
        with open(os.path.join(ROOT_DIR, "no_strand_overlap", "Chr1.gff"), "rb") as f_in:
            with gzip.open(cls.gff_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir, ignore_errors=True)

    def test_create_db_from_gzipped_gtf(self):
        """gffutils.create_db handles a .gtf.gz file correctly."""
        db_path = os.path.join(self.temp_dir, "test_gtf.db")
        gffutils.create_db(self.gtf_gz, db_path, force=True, merge_strategy="create_unique")
        self.assertTrue(os.path.isfile(db_path))
        self.assertGreater(os.path.getsize(db_path), 0)
        db = FeatureDB(db_path)
        # case3.no_parent.gtf has only exons; gffutils infers a gene from them.
        gene = db["ENSMUSG00000112145"]
        self.assertEqual(gene.featuretype, "gene")

    def test_create_db_from_gzipped_gff(self):
        """gffutils.create_db handles a .gff.gz file correctly."""
        db_path = os.path.join(self.temp_dir, "test_gff.db")
        gffutils.create_db(self.gff_gz, db_path, force=True, merge_strategy="create_unique")
        self.assertTrue(os.path.isfile(db_path))
        self.assertGreater(os.path.getsize(db_path), 0)
        db = FeatureDB(db_path)
        # Chr1.gff has an explicit protein_coding_gene feature with ID=C4B63_1g127.
        gene = db["C4B63_1g127"]
        self.assertEqual(gene.featuretype, "protein_coding_gene")


if __name__ == '__main__':
    unittest.main()