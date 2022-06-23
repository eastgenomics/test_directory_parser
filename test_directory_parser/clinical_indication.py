import regex

from test_directory_parser import utils


class ClinicalIndication:
    def __init__(self, r_code, name, panels, test_method):
        self.r_code = r_code
        self.name = name
        self.panels = panels
        self.test_method = test_method
        self.clean_panels()

    def clean_panels(self):
        # stupid weird dash that needs replacing
        self.panels = self.panels.replace("–", "-")
        panels_comma = [p.strip() for p in self.panels.split(",")]
        panels_semicolon = [p.strip() for p in self.panels.split(";")]

        if panels_comma == panels_semicolon:
            # regex to identify panelapp panels
            if utils.match_target(
                r"[A-Za-z0-9-()\ ]*\([0-9&\ ]+\)", panels_comma[0]
            ):
                self.panels = panels_comma
                return

            # regex to identify gene symbol
            if utils.match_target(r"[A-Z]+[A-Z0-9]+", panels_comma[0]):
                self.panels = panels_comma
                return

            # regex to identify the rest
            if utils.match_target(r"[A-Za-z\ ]", panels_comma[0]):
                print(panels_comma)
                self.panels = panels_comma
                return

            print(panels_comma)

        else:
            if len(panels_comma) == 1:
                # try and rescue some panelapp panels
                if utils.match_target(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_comma[0]
                ):
                    self.panels = panels_comma
                    return
                else:
                    # assume that we have lists of genes using semicolon
                    pass
                    # print("assume lists of gene semicolon", panels_comma)

            elif len(panels_comma) >= 2:
                # check if the list only contains genes
                # check that the first element has a gene structure i.e. it's
                # not something weird
                if utils.match_target(r"[A-Z]+[A-Z0-9]+", panels_comma[0]):
                    # gene
                    for panel in panels_comma:
                        if utils.match_target(r"[A-Z]+[A-Z0-9]+", panel):
                            pass
                else:
                    # some panelapp panels have commas
                    self.panels = ", ".join(panels_comma)
                    return

            if len(panels_semicolon) == 1:
                # try and rescue some panelapp panels
                if utils.match_target(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_semicolon[0]
                ):
                    self.panels = panels_semicolon
                    return
                else:
                    # assume that we have lists of genes not using comma
                    # print("assume lists of gene comma", panels_semicolon)
                    pass

            elif len(panels_semicolon) >= 2:
                pass
