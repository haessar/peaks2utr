import os
import os.path
import unittest
from unittest.mock import MagicMock, patch

import gffutils

from peaks2utr import prepare_argparser
from peaks2utr.validation import matching_chr

TEST_DIR = os.path.dirname(__file__)


class TestValidation(unittest.TestCase):
    def setUp(self):
        argparser = prepare_argparser()
        self.args = argparser.parse_args(["Chr1.gtf", ""])
        self.db_path = os.path.join(TEST_DIR, "Chr1.db")
        gffutils.create_db(os.path.join(TEST_DIR, self.args.GFF_IN), self.db_path, force=True)

    def tearDown(self):
        os.remove(self.db_path)

    def test_matching_chr(self):
        mock_af = MagicMock()
        mock_af.fetch.return_value = object
        with patch("peaks2utr.validation.index_bam_file") as mock_index:
            with patch("pysam.AlignmentFile", return_value=mock_af):
                self.assertTrue(matching_chr(self.db_path, self.args))
                mock_af.fetch.side_effect = ValueError()
                self.assertFalse(matching_chr(self.db_path, self.args))
                self.assertEqual(mock_index.call_count, 2)


if __name__ == '__main__':
    unittest.main()
