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

    def test_parse_config(self):
        """ Test for parse_config function """

        path_to_test_data = PATH_TO_TEST_FOLDER / "test.json"
        test_output = parse_config(path_to_test_data)
        expected_output = {"test_field": "test_value"}
        self.assertEqual(test_output, expected_output)

    @patch("test_directory_parser.rare_disease.pd.read_excel")
    def test_parse_rare_disease_td(self, mock_td_excel):
        """ Test for parse_rare_disease_td """

        test_config = {
            "name": "Test_config",
            "sheet_of_interest": "R&ID indications",
            "clinical_indication_column_code": "Test ID",
            "clinical_indication_column_name": "Clinical Indication",
            "panel_column": "Target/Genes",
            "test_method_column": "Test Method",
            "changes_column": "Changes",
            "specialty_column": "Specialist test group",
            "header_index": 1,
            "specialism_of_interest": ["Core", "Endocrinology", "Neurology"],
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

        mock_td_excel.return_value = {
            "R&ID indications": pd.DataFrame(
                {
                    "Clinical indication ID": ["R100"],
                    "Test ID": ["R100.1"],
                    "Clinical Indication": ["CI1"],
                    "Target/Genes": ["Panel1"],
                    "Test Method": ["WES"],
                    "Commissioning category": ["Category1"],
                    "Specialist test group": ["Core"],
                    "Changes since April 2023 publication": ["No change"]
                }
            )
        }

        expected_output = pd.DataFrame(
            {
                "Test ID": ["R100.1"],
                "Clinical Indication": ["CI1"],
                "Target/Genes": ["Panel1"],
                "Test Method": ["WES"],
                "Specialist test group": ["Core"],
                "Changes since April 2023 publication": ["No change"]
            }
        )

        test_output, test_change_column = parse_rare_disease_td("", test_config)

        with self.subTest():
            self.assertTrue(test_output.equals(expected_output))
        
        with self.subTest():
            self.assertEqual(
                test_change_column, "Changes since April 2023 publication"
            )
