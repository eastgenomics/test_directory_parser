import pandas as pd
import regex

from .test_directory_parser import utils


class ClinicalIndication:
    def __init__(
        self, r_code: str, name: str, target: str, test_method: str,
        change: str, hgnc_dump: pd.DataFrame
    ):
        """ Setup the clinical indication object

        Args:
            r_code (str): R code associated with the clinical indication
            name (str): Name of the clinical indication
            target (str): Target of the clinical indication
            test_method (str): Test method of the clinical indication
            change (str): Changes described in the test directory excel i.e.
            no changes, addition of, removal of ...
            hgnc_dump (pd.DataFrame): Dataframe containing the HGNC data
        """

        self.r_code = r_code
        # R424 has a space at the end of its name
        self.name = name.replace("–", "-").strip()

        if "gene" in test_method:
            self.gemini_name = f"{self.r_code}_{self.name}_G"
        else:
            self.gemini_name = f"{self.r_code}_{self.name}_P"

        self.original_targets = target
        self.test_method = test_method
        self.change = change
        self.clean_target(hgnc_dump)

    def clean_target(self, hgnc_dump: pd.DataFrame):
        """ Attempt to clean up the targets in the excel file

        Args:
            hgnc_dump (pandas.Dataframe): Dataframe of hgnc data
        """

        self.panels = []
        self.genes = []

        potential_panel_targets = regex.findall(
            r"\([0-9&\ ]+\)", self.original_targets
        )
        potential_gene_targets = regex.findall(
            r"[A-Z]+[A-Z0-9\-]+", self.original_targets
        )

        # regex to identify panelapp panels
        if potential_panel_targets:
            for potential_panel in potential_panel_targets:
                cleaned_panelapp_id = potential_panel.replace(
                    "(", "").replace(")", "")
                self.panels.append(cleaned_panelapp_id)

        # regex to identify gene symbol
        elif potential_gene_targets:
            for potential_gene in potential_gene_targets:
                hgnc_id_data = utils.find_hgnc_id(potential_gene, hgnc_dump)

                if hgnc_id_data["HGNC ID"]:
                    self.genes.append(hgnc_id_data["HGNC ID"])
                else:
                    self.genes.append(None)
