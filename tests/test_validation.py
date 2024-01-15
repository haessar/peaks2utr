import os
import os.path
import unittest
from unittest.mock import MagicMock, patch

import gffutils

from peaks2utr.validation import matching_chr

TEST_DIR = os.path.dirname(__file__)


class TestValidation(unittest.TestCase):
    def setUp(self):
        self.db_path = os.path.join(TEST_DIR, "Chr1.db")
        gffutils.create_db(os.path.join(TEST_DIR, "Chr1.gtf"), self.db_path, force=True)

    def tearDown(self):
        os.remove(os.path.join(TEST_DIR, "Chr1.db"))

    def test_matching_chr(self):
        mock_af = MagicMock()
        mock_af.fetch.return_value = object
        with patch("pysam.AlignmentFile", return_value=mock_af):
            self.assertTrue(matching_chr(self.db_path, ""))
            mock_af.fetch.side_effect = ValueError()
            self.assertFalse(matching_chr(self.db_path, ""))


if __name__ == '__main__':
    unittest.main()
