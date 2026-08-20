import asyncio
import os.path
import unittest
from unittest.mock import patch

from peaks2utr import prepare_argparser
from peaks2utr.preprocess import create_db


class TestCreateDBOptions(unittest.TestCase):
    def test_infer_options_default_to_false(self):
        args = prepare_argparser().parse_args(["annotations.gtf", "reads.bam"])

        self.assertFalse(args.disable_infer_genes)
        self.assertFalse(args.disable_infer_transcripts)

    def test_infer_options_are_parsed(self):
        args = prepare_argparser().parse_args([
            "annotations.gtf",
            "reads.bam",
            "--disable-infer-genes",
            "--disable-infer-transcripts",
        ])

        self.assertTrue(args.disable_infer_genes)
        self.assertTrue(args.disable_infer_transcripts)

    def test_create_db_passes_infer_options_to_gffutils(self):
        db_path = os.path.join("cache", "annotations.db")

        with patch("peaks2utr.preprocess.cached", return_value=db_path):
            with patch("peaks2utr.preprocess.os.path.isfile", return_value=False):
                with patch("peaks2utr.preprocess.gffutils.create_db") as mock_create_db:
                    result = asyncio.run(create_db(
                        "annotations.gtf",
                        disable_infer_genes=True,
                        disable_infer_transcripts=True,
                    ))

        self.assertEqual(result, db_path)
        mock_create_db.assert_called_once_with(
            "annotations.gtf",
            db_path,
            force=True,
            verbose=True,
            merge_strategy="create_unique",
            disable_infer_genes=True,
            disable_infer_transcripts=True,
        )


if __name__ == '__main__':
    unittest.main()
