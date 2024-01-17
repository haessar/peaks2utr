import os
import os.path
import unittest
from unittest.mock import MagicMock, patch

from peaks2utr import prepare_argparser
from peaks2utr.exceptions import PysamError
from peaks2utr.validation import matching_chr, valid_bam

TEST_DIR = os.path.dirname(__file__)


class TestValidation(unittest.TestCase):
    def setUp(self):
        argparser = prepare_argparser()
        self.args = argparser.parse_args([os.path.join(TEST_DIR, "Chr1.gtf"), ""])

    def test_matching_chr(self):
        mock_af = MagicMock()
        mock_af.fetch.return_value = object
        with patch("peaks2utr.validation.index_bam_file") as mock_index:
            with patch("pysam.AlignmentFile", return_value=mock_af):
                self.assertTrue(matching_chr(self.args))
                mock_af.fetch.side_effect = ValueError()
                self.assertFalse(matching_chr(self.args))
                self.assertEqual(mock_index.call_count, 2)

    def test_valid_bam(self):
        mock_af = MagicMock()
        with patch("builtins.open"):
            with patch("pysam.AlignmentFile", return_value=mock_af):
                self.assertTrue(valid_bam(self.args))
            with patch("pysam.AlignmentFile", side_effect=ValueError):
                with self.assertRaises(PysamError):
                    valid_bam(self.args)


if __name__ == '__main__':
    unittest.main()
