import json
from pathlib import Path

import pandas as pd

from .test_directory_parser import clinical_indication
from .test_directory_parser import rare_disease
from .test_directory_parser import utils


class TestDirectory:
    def __init__(
        self, test_directory_path: str, config_path: str, td_type: str,
        hgnc_dump: pd.DataFrame
    ):
        """ Setup the test directory object with its clinical indications

        Args:
            test_directory_path (str): Path to the test directory
            config_path (str): Path to the config file for parsing the test
            directory
            td_type (str): Type of the test directory i.e. rare disease, cancer
            hgnc_dump (pd.DataFrame): Dataframe containing the HGNC data for
            symbols conversion
        """

        config_data = rare_disease.parse_config(config_path)
        sheet, change_column = rare_disease.parse_rare_disease_td(
            test_directory_path, config_data
        )

        # convert the dataframe to a dictionary
        self.data = sheet.to_dict()
        self.td = Path(test_directory_path).name
        self.td_type = td_type
        self.change_column = change_column
        self.config = config_data
        self.all_clinical_indications = []
        self.ngs_clinical_indications = []
        self.filtered_clinical_indications = []
        self.hgnc_dump = hgnc_dump

    def setup_clinical_indications(self):
        """ Go through the test directory and create clinical indications """

        r_codes = self.data[
            self.config["clinical_indication_column_code"]
        ]
        clinical_indications = self.data[
            self.config["clinical_indication_column_name"]
        ]
        panels = self.data[
            self.config["panel_column"]
        ]
        test_methods = self.data[
            self.config["test_method_column"]
        ]
        changes = self.data[self.change_column]

        for index, ci in clinical_indications.items():
            r_code = r_codes[index]
            panel = panels[index]
            test_method = test_methods[index]
            change = changes[index]

            ci = clinical_indication.ClinicalIndication(
                r_code, ci, panel, test_method, change, self.hgnc_dump
            )

            # handled clinical indications by the lab and that will be stored
            # in panel palace
            if test_method.strip() in self.config["ngs_test_methods"]:
                self.ngs_clinical_indications.append(ci)

            self.all_clinical_indications.append(ci)

    def filter_clinical_indications(self, internal_td: pd.DataFrame):
        """ Filter the gathered clinical indications using the internal test
        directory

        Args:
            internal_td (pd.DataFrame): Dataframe containing the internal test
            directory data
        """

        for clinical_indication in self.ngs_clinical_indications:
            filtered_df = internal_td[
                internal_td["Test ID"] == clinical_indication.r_code
            ]

            if filtered_df.shape[0] != 0:
                self.filtered_clinical_indications.append(clinical_indication)

        assert (
            len(self.filtered_clinical_indications) == internal_td.shape[0]
        ), (
            "Did not found every test from the internal test directory in the "
            "new test directory"
        )

    def output_json(self, output: str):
        """ Output the content of the test directory object as a JSON file

        Args:
            output (str): Path to the output
        """

        print("\nOutputting json file..\n")
        td = self.td
        source = self.config["name"]
        date = utils.get_date()

        indications = []

        if self.filtered_clinical_indications:
            clinical_indications = self.filtered_clinical_indications
        else:
            clinical_indications = self.ngs_clinical_indications

        for ci in clinical_indications:
            # we only output clinical indications that the lab will handle
            indication = {
                "name": ci.name, "code": ci.r_code,
                "gemini_name": ci.gemini_name, "test_method": ci.test_method,
                "panels": ci.panels + ci.genes,
                "original_targets": ci.original_targets, "changes": ci.change
            }

            if None in indication["panels"]:
                print((
                    f"Check {ci.r_code} for why the target is or contains "
                    f"None: {ci.original_targets}"
                ))

            indications.append(indication)

        data = {
            "td_source": td, "config_source": source, "date": date,
            "indications": indications
        }

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
