import unittest
from unittest.mock import patch

from peaks2utr import main


class TestEntryPoint(unittest.TestCase):
    def test_main(self):
        with patch('asyncio.run') as mock_run:
            with patch('peaks2utr.prepare_argparser') as mock_argparser:
                main()
                self.assertEqual(mock_run.call_count, 1)
                self.assertEqual(mock_argparser.call_count, 1)


if __name__ == '__main__':
    unittest.main()
