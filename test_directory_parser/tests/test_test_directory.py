import json
from pathlib import Path
import unittest
from unittest.mock import patch

import pandas as pd

from test_directory_parser.clinical_indication import ClinicalIndication
from test_directory_parser.test_directory import TestDirectory
from test_directory_parser.utils import get_date

PATH_TO_TEST_FOLDER = Path(".") / "test_directory_parser" / "tests" / "test_data"

class TestTestDirectory(unittest.TestCase):
    """ Test the methods of the test directory object """

    @classmethod
    @patch("test_directory_parser.rare_disease.parse_rare_disease_td")
    @patch("test_directory_parser.rare_disease.parse_config")
    def setUpClass(cls, mock_config_data, mock_td_data):
        """ Test the Test directory object """

        mock_config_data.return_value = {
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

        mock_td_data.return_value = (
            pd.DataFrame(
                {
                    "Test ID": ["R100.1", "R200.1"],
                    "Clinical Indication": ["CI1", "CI2"],
                    "Target/Genes": ["Panel1", "Panel2"],
                    "Test Method": ["WES", "Not WES"],
                    "Specialist test group": ["Core", "Core"],
                    "Changes since April 2023 publication": ["No change", "No change"]
                }
            ),
            "Changes since April 2023 publication"
        )

        td = TestDirectory("path/to/test_test_directory", "", "rare_disease", "")
        td.setup_clinical_indications()
        cls.td = td

    def tearDown(self) -> None:
        # reset the filtered_clinical_indications attributes for further testing 
        self.td.filtered_clinical_indications = []

    def test_setup_clinical_indications(self):
        """ Test for setting up clinical indications i.e. adding them to the
        NGS clinical indications or not
        
        Test setup:
        - TD contains R100.1 and R200.1
        - R100.1 is WES
        - R200.1 is not WES

        Expectation:
        - R100.1 in ngs_clinical_indications
        - R100.1 and R200.1 in all_clinical_indications 
        """

        # the clinical indications variables are lists containing clinical
        # indication objects and I don't have ways to compare those objects
        # instead, I check the number of objects in each lists
        with self.subTest():
            self.assertEqual(len(self.td.ngs_clinical_indications), 1)

        with self.subTest():
            self.assertEqual(len(self.td.all_clinical_indications), 2)

    def test_filter_clinical_indications_present(self):
        """ Test filter clinical indications method with internal TD having a
        clinical indication 

        Test setup:
        - TD contains R100.1 and R200.1
        - Internal TD contains R100.1

        All the tests in the internal TD are present in the "normal" TD, check
        if they are added correctly in the filtered_clinical_indications
        attribute
        """

        test_internal_td = pd.DataFrame(
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

        self.td.filter_clinical_indications(test_internal_td)
        self.assertEqual(len(self.td.filtered_clinical_indications), 1)

    def test_filter_clinical_indications_absent(self):
        """ Test filter clinical indications method with internal TD missing a
        clinical indication
        
        Test setup:
        - TD contains R100.1 and R200.1
        - Internal TD contains R100.1 and R100.2

        The code expects the tests in the internal TD to be all present in the
        "normal" TD triggering the error
        """

        test_internal_td = pd.DataFrame(
            {
                "Clinical indication ID": ["R100", "R100"],
                "Test ID": ["R100.1", "R100.2"],
                "Clinical Indication": ["CI1", "CI2"],
                "Target/Genes": ["Panel1", "Panel2"],
                "Test Method": ["WES", "WES"],
                "Commissioning category": ["Category1", "Category2"],
                "Specialist test group": ["Core", "Core"],
                "Changes since April 2023 publication": ["No change", "No change"]
            }
        )

        with self.assertRaises(AssertionError):
            self.td.filter_clinical_indications(test_internal_td)

    def test_output_json_no_internal_td(self):
        """ Test output json method
        
        Test setup:
        - TD object contains only R100.1
        - No internal TD filtered

        Expectation:
        - Output JSON should contain data only for R100.1
        """

        expected_output = {
            "td_source": "test_test_directory",
            "config_source": "Test_config",
            "date": get_date(),
            "indications": [
                {
                "name": "CI1",
                "code": "R100.1",
                "gemini_name": "R100.1_CI1_P",
                "test_method": "WES",
                "panels": [],
                "original_targets": "Panel1",
                "changes": "No change"
                }
            ]
        }

        self.td.output_json(PATH_TO_TEST_FOLDER / "test_output.json")

        with open(PATH_TO_TEST_FOLDER / "test_output.json") as f:
            self.assertEqual(json.loads(f.read()), expected_output)
            

    @patch("test_directory_parser.utils.find_hgnc_id")
    def test_output_json_with_internal_td(self, mock_hgnc_id):
        """ Test output json method
        
        Test setup:
        - TD object contains only R100.1
        - No internal TD filtered

        Expectation:
        - Output JSON should contain data only for R100.1
        """

        # mock in a return value so that the gene symbol is not identified
        mock_hgnc_id.return_value = pd.Series(
            {
                "HGNC ID": None
            }
        )

        # create dummy clinical indication with unknown gene symbol
        test_ci = ClinicalIndication(
            "R5.1", "Unknown gene CI", "BLARG", "Single gene panel", "", ""
        )

        # bypass filtering of clinical indications by assigning to the
        # attribute directly 
        self.td.filtered_clinical_indications = [test_ci]

        expected_output = {
            "td_source": "test_test_directory",
            "config_source": "Test_config",
            "date": get_date(),
            "indications": [
                {
                "name": "Unknown gene CI",
                "code": "R5.1",
                "gemini_name": "R5.1_Unknown gene CI_G",
                "test_method": "Single gene panel",
                "panels": [None],
                "original_targets": "BLARG",
                "changes": ""
                }
            ]
        }

        self.td.output_json(PATH_TO_TEST_FOLDER / "test_output.json")

        with open(PATH_TO_TEST_FOLDER / "test_output.json") as f:
            self.assertEqual(json.loads(f.read()), expected_output)
            
