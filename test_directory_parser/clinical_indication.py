import regex

from test_directory_parser import utils


class ClinicalIndication:
    def __init__(self, r_code, name, panels, test_method, hgnc_dump):
        self.r_code = r_code
        self.name = name
        self.original_targets = panels
        self.panels = None
        self.test_method = test_method
        self.clean_panels(hgnc_dump)

    def clean_panels(self, hgnc_dump):
        """ Attempt to clean up the targets in the excel file

        Args:
            hgnc_dump (pandas.Dataframe): Dataframe of hgnc data
        """

        # stupid weird dash that needs replacing
        panels = self.original_targets.replace("–", "-")
        panels_comma = [p.strip() for p in panels.split(",")]
        panels_semicolon = [p.strip() for p in panels.split(";")]

        if panels_comma == panels_semicolon:
            # regex to identify panelapp panels
            if regex.match(
                r"[A-Za-z0-9-()\ ]*\([0-9&\ ]+\)", panels_comma[0]
            ):
                self.panels = panels_comma
                return

            # regex to identify gene symbol
            if regex.match(r"[A-Z]+[A-Z0-9]+", panels_comma[0]):
                hgnc_id = utils.find_hgnc_id(panels_comma[0], hgnc_dump)

                if hgnc_id:
                    self.panels = [hgnc_id]

                return

            # regex to identify the rest
            if regex.match(r"[A-Za-z\ ]", panels_comma[0]):
                return

        else:
            if len(panels_comma) == 1:
                # try and rescue some panelapp panels
                if regex.match(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_comma[0]
                ):
                    self.panels = panels_comma
                    return
                else:
                    # assume that we have lists of genes using semicolon
                    pass
                    # print("assume lists of gene semicolon", panels_comma)

            elif len(panels_comma) >= 2:
                cleaned_panels = utils.handle_list_panels(
                    panels_comma, hgnc_dump, self.r_code
                )

                if cleaned_panels:
                    self.panels = cleaned_panels
                return

            if len(panels_semicolon) == 1:
                # try and rescue some panelapp panels
                if regex.match(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_semicolon[0]
                ):
                    self.panels = panels_semicolon
                    return
                else:
                    # assume that we have lists of genes not using comma
                    # print("assume lists of gene comma", panels_semicolon)
                    pass

            elif len(panels_semicolon) >= 2:
                cleaned_panels = utils.handle_list_panels(
                    panels_semicolon, hgnc_dump, self.r_code
                )

                if cleaned_panels:
                    self.panels = cleaned_panels
                return
