from test_directory_parser import clinical_indication


class TestDirectory:
    def __init__(self, dataframe, config, td_type):
        self.data = dataframe.to_dict()
        self.type = td_type
        self.config = config
        self.clinical_indications = []

    def setup_clinical_indications(self):
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

            self.clinical_indications.append(
                clinical_indication.ClinicalIndication(
                    r_code, ci, panel, test_method
                )
            )
