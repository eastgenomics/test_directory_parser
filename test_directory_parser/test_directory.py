import json
from pathlib import Path

from test_directory_parser import clinical_indication
from test_directory_parser import rare_disease
from test_directory_parser import utils


class TestDirectory:
    def __init__(
            self, test_directory_path, config_path, td_type, hgnc_dump
        ):
        config_data = rare_disease.parse_config(config_path)
        sheet, change_column = rare_disease.parse_rare_disease_td(
            test_directory_path, config_data
        )
        self.data = sheet.to_dict()
        self.td = Path(test_directory_path).name
        self.type = td_type
        self.change_column = change_column
        self.config = config_data
        self.all_clinical_indications = []
        self.ngs_clinical_indications = []
        self.hgnc_dump = hgnc_dump

    def setup_clinical_indications(self):
        """ Setup the clinical indications """

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

    def output_json(self, output):
        print("\nOutputting json file..\n")
        td = self.td
        source = self.config["name"]
        date = utils.get_date()

        indications = []

        for ci in self.ngs_clinical_indications:
            # we only output clinical indications that the lab will handle
            indication = {
                "name": ci.name, "code": ci.r_code,
                "gemini_name": ci.gemini_name, "test_method": ci.test_method,
                "panels": ci.panels, "original_targets": ci.original_targets,
                "changes": ci.change
            }

            if ci.panels is None or None in ci.panels:
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
