import regex


class ClinicalIndication:
    def __init__(self, r_code, name, panels, test_method):
        self.r_code = r_code
        self.name = name
        self.panels = panels
        self.test_method = test_method
        self.clean_panels()

    def clean_panels(self):
        self.panels = self.panels.replace("–", "-")
        panels_comma = [p.strip() for p in self.panels.split(",")]
        panels_semicolon = [p.strip() for p in self.panels.split(";")]

        if panels_comma == panels_semicolon:
            # regex to identify panelapp panels
            match = regex.match(r"[A-Za-z0-9-()\ ]*\([0-9]+\)", panels_comma[0])

            # the panel is panelapp panel
            if match:
                self.panels = panels_comma[0]
                return

            # regex to identify gene symbol
            match = regex.match(r"[A-Z]+[A-Z0-9]+", panels_comma[0])

            if match:
                self.panels = panels_comma[0]
                return

            # regex to identify the rest
            match = regex.match(r"[A-Za-z\ ]", panels_comma[0])

            if match:
                self.panels = panels_comma[0]
                return

        else:
            if len(panels_comma) == 1:
                # try and rescue some panelapp panels
                match = regex.match(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_comma[0]
                )

                if match:
                    self.panels = panels_comma[0]
                    return

            elif len(panels_comma) >= 2:
                pass

            if len(panels_semicolon) == 1:
                # try and rescue some panelapp panels
                match = regex.match(
                    r"[A-Za-z0-9-()\ ,]*\([0-9]+\)", panels_semicolon[0]
                )

                if match:
                    self.panels = panels_semicolon[0]
                    return

            elif len(panels_semicolon) >= 2:
                pass
