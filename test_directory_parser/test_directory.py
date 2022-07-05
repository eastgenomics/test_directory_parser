from test_directory_parser import clinical_indication


class TestDirectory:
    def __init__(self, dataframe, config, td_type, hgnc_dump):
        self.data = dataframe.to_dict()
        self.type = td_type
        self.config = config
        self.clinical_indications = []
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

            self.clinical_indications.append(
                clinical_indication.ClinicalIndication(
                    r_code, ci, panel, test_method, self.hgnc_dump
                )
            )

    def output(self):
        """ Output for every clinical indication in the test methods of
        interest: r_code, name, test_method, original targets and cleaned
        targets
        """

        with open("cleaned_test_directory.tsv", "w") as f:
            for ci in self.clinical_indications:
                if ci.test_method in self.config["ngs_test_methods"]:
                    if ci.panels is None:
                        f.write(
                            f"{ci.r_code}\t{ci.name}\t{ci.test_method}\t"
                            f"{ci.original_targets}\t{ci.panels}\n"
                        )
                    else:
                        f.write(
                            f"{ci.r_code}\t{ci.name}\t{ci.test_method}\t"
                            f"{ci.original_targets}\t{'|'.join(ci.panels)}\n"
                        )
