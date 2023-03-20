import json

from test_directory_parser import clinical_indication
from test_directory_parser import utils


class TestDirectory:
    def __init__(self, dataframe, config, td_type, hgnc_dump):
        self.data = dataframe.to_dict()
        self.type = td_type
        self.config = config
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

        for index, ci in clinical_indications.items():
            r_code = r_codes[index]
            panel = panels[index]
            test_method = test_methods[index]

            ci = clinical_indication.ClinicalIndication(
                r_code, ci, panel, test_method, self.hgnc_dump
            )

            # handled clinical indications by the lab and that will be stored
            # in panel palace
            if test_method.strip() in self.config["ngs_test_methods"]:
                self.ngs_clinical_indications.append(ci)

            self.all_clinical_indications.append(ci)

    def output_json(self, output):
        source = self.config["name"]
        date = utils.get_date()

        # we only output clinical indications that the lab will handle
        indications = [
            {
                "name": ci.name, "code": ci.r_code,
                "gemini_name": ci.gemini_name, "test_method": ci.test_method,
                "panels": ci.panels, "original_targets": ci.original_targets
            }
            for ci in self.ngs_clinical_indications
        ]

        data = {
            "source": source, "date": date, "indications": indications
        }

        with open(output, "w") as f:
            json.dump(data, f, indent=2)
