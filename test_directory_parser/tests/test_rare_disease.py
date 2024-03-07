from pathlib import Path
import unittest
from unittest.mock import patch

import pandas as pd

from test_directory_parser.rare_disease import (
    parse_config, parse_rare_disease_td
)

PATH_TO_TEST_FOLDER = Path(".") / "test_directory_parser" / "tests" / "test_data"

class TestRareDisease(unittest.TestCase):
    """ Suite of tests for the rare_disease.py script """

    def setUp(self) -> None:
        self.test_config = {
            "name": "Test config",
            "sheet_of_interest": "Sheet1",
            "clinical_indication_column_code": "Test ID",
            "clinical_indication_column_name": "Clinical Indication",
            "panel_column": "Target/Genes",
            "test_method_column": "Test Method",
            "changes_column": "Changes",
            "ngs_column": "NGS Technology",
            "header_index": 0,
            "ngs_type": ["WES", "CEN"],
            "ngs_test_methods": [
                "Medium panel", "Single gene sequencing <=10 amplicons",
                "Single gene sequencing <10 amplicons",
                "Single gene sequencing >=10 amplicons",
                "Single gene testing (<10 amplicons)", "small panel", "Small panel",
                "WES or Large panel", "WES or Large Panel", "WES or Medium Panel",
                "WES or Large penel", "WES or Medium panel", "WES or Small Panel",
                "WGS", "WES", "Exon level CNV detection by MLPA or equivalent"
            ]
        }

    def test_parse_config(self):
        """ Test for parse_config function """

        path_to_test_data = PATH_TO_TEST_FOLDER / "test.json"
        test_output = parse_config(path_to_test_data)
        expected_output = {"test_field": "test_value"}
        self.assertEqual(test_output, expected_output)

    @patch("test_directory_parser.rare_disease.pd.read_excel")
    def test_parse_rare_disease_td_correctly(self, mock_td_excel):
        """ Test for parse_rare_disease_td 
        
        Setup test:
        - Normal test directory

        Expectations:
        - Found change column
        - Correct parsing of td_excel data
        """

        mock_td_excel.return_value = pd.DataFrame(
            {
                "Clinical indication ID": ["R100", "R300"],
                "Test ID": ["R100.1", "R300.1"],
                "Clinical Indication": ["CI1", "CI2"],
                "Target/Genes": ["Panel1", "Panel2"],
                "Test Method": ["WES", "WES"],
                "Commissioning category": ["Category1", "Category2"],
                "NGS Technology": ["WES", "Not WES"],
                "Changes since April 2023 publication": ["No change", "No change"]
            }
        )

        expected_output = pd.DataFrame(
            {
                "Test ID": ["R100.1"],
                "Clinical Indication": ["CI1"],
                "Target/Genes": ["Panel1"],
                "Test Method": ["WES"],
                "NGS Technology": ["WES"],
                "Changes since April 2023 publication": ["No change"]
            }
        )

        test_output, test_change_column = parse_rare_disease_td(
            "", self.test_config
        )

        with self.subTest():
            self.assertTrue(test_output.equals(expected_output))
        
        with self.subTest():
            self.assertEqual(
                test_change_column, "Changes since April 2023 publication"
            )

    @patch("test_directory_parser.rare_disease.pd.read_excel")
    def test_parse_rare_disease_td_too_many_change_columns(
        self, mock_td_excel
    ):
        """ Test for parse_rare_disease_td

        Setup test:
        - Test directory with 2 columns with "changes" in their header

        Expectations:
        - Exception raised because 2 >= columns with changes in their header 
        """

        mock_td_excel.return_value = pd.DataFrame(
            {
                "Clinical indication ID": ["R100"],
                "Test ID": ["R100.1"],
                "Clinical Indication": ["CI1"],
                "Target/Genes": ["Panel1"],
                "Test Method": ["WES"],
                "Commissioning category": ["Category1"],
                "Specialist test group": ["Core"],
                "Changes since April 2023 publication": ["No change"],
                "Other changes": ["Kinda change"]
            }
        )

        with self.assertRaises(Exception):
            test_output, test_change_column = parse_rare_disease_td(
                "", self.test_config
            )

    @patch("test_directory_parser.rare_disease.pd.read_excel")
    def test_parse_rare_disease_td_no_change_column(
        self, mock_td_excel
    ):
        """ Test for parse_rare_disease_td

        Setup test:
        - Test directory with 0 columns with "changes" in their header

        Expectations:
        - Exception raised because 0 columns with "changes" in their header 
        """

        mock_td_excel.return_value = pd.DataFrame(
            {
                "Clinical indication ID": ["R100"],
                "Test ID": ["R100.1"],
                "Clinical Indication": ["CI1"],
                "Target/Genes": ["Panel1"],
                "Test Method": ["WES"],
                "Commissioning category": ["Category1"],
                "Specialist test group": ["Core"],
            }
        )

        with self.assertRaises(Exception):
            test_output, test_change_column = parse_rare_disease_td(
                "", self.test_config
            )
