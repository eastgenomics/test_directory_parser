import unittest
from unittest.mock import patch

import pandas as pd

from test_directory_parser.clinical_indication import ClinicalIndication

class TestClinicalIndication(unittest.TestCase):
    """ Class to test the clean_target method of the ClinicalIndication object """
    multiple_genes_side_effect=[
        pd.Series(
            {
                "Gene symbol": "BRCA1",
                "HGNC ID": "HGNC:1100",
                "Previous": None,
                "Alias": None
            }
        ),
        pd.Series(
            {
                "Gene symbol": "BRCA2",
                "HGNC ID": "HGNC:1101",
                "Previous": None,
                "Alias": None
            }
        )
    ]

    def test_clean_target_one_panel(self):
        """ Test the clean_target method for one panel """
        test_clinical_indication = ClinicalIndication(
            "R100.1", "CI1", "Panelapp panel 1 (100)", "WES",
            "Things have changed", ""
        )

        test_output = test_clinical_indication.panels
        expected_output = ["100"]

        self.assertEqual(test_output, expected_output)

    def test_clean_target_multiple_panels(self):
        """ Test the clean_target method for multiple panels """
        test_clinical_indication = ClinicalIndication(
            "R100.1", "CI1", "Panelapp panel 1 (100), Panelapp panel 2 (101)",
            "WES", "Things have changed", ""
        )

        test_output = test_clinical_indication.panels
        expected_output = ["100", "101"]

        self.assertEqual(test_output, expected_output)

    @patch("test_directory_parser.clinical_indication.utils.find_hgnc_id")
    def test_clean_target_one_gene(self, mock_hgnc_data):
        """ Test the clean_target method for one gene """
        mock_hgnc_data.return_value = pd.Series(
            {
                "Gene symbol": "BRCA1",
                "HGNC ID": "HGNC:1100",
                "Previous": None,
                "Alias": None
            }
        )
        test_clinical_indication = ClinicalIndication(
            "R100.1", "CI1", "BRCA1", "WES",
            "Things have changed", ""
        )

        test_output = test_clinical_indication.genes
        expected_output = ["HGNC:1100"]

        self.assertEqual(test_output, expected_output)

    @patch(
        "test_directory_parser.clinical_indication.utils.find_hgnc_id",
        side_effect=multiple_genes_side_effect
    )
    def test_clean_target_multiple_genes(self, mock_hgnc_data):
        """ Test the clean_target method for multiple genes """
        test_clinical_indication = ClinicalIndication(
            "R100.1", "CI1", "BRCA1, BRCA2", "WES",
            "Things have changed", ""
        )

        test_output = test_clinical_indication.genes
        expected_output = ["HGNC:1100", "HGNC:1101"]

        self.assertEqual(test_output, expected_output)
        self.assertTrue(mock_hgnc_data.call_count, 2)

    def test_clean_target_mixin(self):
        """ Test the clean_target method for a mix of panels and genes """
        test_clinical_indication = ClinicalIndication(
            "R100.1", "CI1",
            "Panelapp panel 1 (100), Panelapp panel 2 (101), BRCA1, BRCA2", "WES",
            "Things have changed", ""
        )

        test_gene_output = test_clinical_indication.genes
        expected_gene_output = []
        test_panel_output = test_clinical_indication.panels
        expected_panel_output = ["100", "101"]

        with self.subTest():
            self.assertEqual(test_gene_output, expected_gene_output)

        with self.subTest():
            self.assertEqual(test_panel_output, expected_panel_output)
