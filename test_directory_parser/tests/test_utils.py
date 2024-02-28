import datetime
import numpy as np
import pandas as pd
from pathlib import Path
import unittest
from unittest.mock import patch

from test_directory_parser.utils import (
    get_date, parse_tsv, parse_lab_excel, find_hgnc_id
)

PATH_TO_TEST_FOLDER = Path(".") / "test_directory_parser" / "tests" / "test_data"

class TestUtils(unittest.TestCase):
    def test_get_date(self):
        expected_date = datetime.datetime.now().strftime("%y%m%d")
        test_date = get_date()
        self.assertEqual(expected_date, test_date)

    def test_parse_tsv(self):
        path_to_test_data = PATH_TO_TEST_FOLDER / "test.tsv"
        test_output = parse_tsv(path_to_test_data)
        expected_output = pd.DataFrame(
            {
                "column1": ["test_column1"], "column2": ["test_column2"],
                "column3": ["test_column3"], "column4": [np.nan]
            },
        )
        self.assertTrue(test_output.equals(expected_output))

    @patch("test_directory_parser.utils.pd.read_excel")
    def test_parse_lab_excel(self, mock_excel):
        mock_excel.return_value = pd.DataFrame(
            {
                "Test": ["Test1", "Test2", "Test3", "Test4"],
                "NGS Technology": ["CEN", "WES", "MLPA", "MYE"]
            }
        )
        test_output = parse_lab_excel("")
        expected_output = pd.DataFrame(
            {
                "Test": ["Test1", "Test2"],
                "NGS Technology": ["CEN", "WES"]
            }
        )
        self.assertTrue(test_output.equals(expected_output))

    def test_find_hgnc_id(self):
        test_hgnc_dump = pd.DataFrame(
            {
                "HGNC ID": [
                    "HGNC:1100", "HGNC:28470", "HGNC:1550", "HGNC:1601",
                    "HGNC:11577", "HGNC:24042", "HGNC:13030"
                ],
                "Approved symbol": [
                    "BRCA1", "BRCA1P1", "CBS", "RYR1", "TAFAZZIN", "WWTR1",
                    "ZBTB18"
                ],
                "Status": [
                    "Approved", "Approved", "Approved", "Approved", "Approved",
                    "Approved", "Approved", 
                ],
                "Previous symbols": [
                    "", "", "", "MHS, MHS1, CCO", "CMD3A, EFE2, EFE, TAZ", "",
                    "ZNF238"
                ],
                "Alias symbols": [
                    "RNF53, BRCC1, PPP1R53, FANCS",
                    "LBRCA1, PsiBRCA1, pseudo-BRCA1", "HIP4", "RYR, PPP1R137",
                    "BTHS, G4.5, TAZ1", "TAZ, DKFZp586I1419",
                    "C2H2-171, TAZ-1, RP58", 
                ]
            }
        )

        test_inputs = {
            "BRCA1": pd.Series(
                {
                    "Gene symbol": "BRCA1",
                    "HGNC ID": "HGNC:1100",
                    "Previous": None,
                    "Alias": None
                }
            ).T,
            "TAZ": pd.Series(
                {
                    "Gene symbol": "TAZ",
                    "HGNC ID": None,
                    "Previous": True,
                    "Alias": True
                }
            ),
            "HIP4": pd.Series(
                {
                    "Gene symbol": "HIP4",
                    "HGNC ID": "HGNC:1550",
                    "Previous": False,
                    "Alias": True
                }
            ),
            "Unknown": pd.Series(
                {
                    "Gene symbol": "Unknown",
                    "HGNC ID": None,
                    "Previous": None,
                    "Alias": None
                }
            ),
            "CCO": pd.Series(
                {
                    "Gene symbol": "CCO",
                    "HGNC ID": "HGNC:1601",
                    "Previous": True,
                    "Alias": False
                }
            ),
        }

        for test_input, expected_output in test_inputs.items():
            test_output = find_hgnc_id(test_input, test_hgnc_dump)

            with self.subTest():
                self.assertTrue(test_output.equals(expected_output))


if __name__ == '__main__':
    unittest.main()
