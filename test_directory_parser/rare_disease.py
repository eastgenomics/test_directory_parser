import json

import pandas as pd
from termcolor import colored, cprint


def parse_config(config):
    config_content = open(config)
    config_data = json.load(config_content)
    config_content.close()
    return config_data


def parse_rare_disease_td(test_directory, config):
    xls = pd.read_excel(
        test_directory, sheet_name=None, header=config["header_index"]
    )

    for sheet in xls:
        if sheet == config["sheet_of_interest"]:
            data = xls[sheet].loc(axis=1)[
                config["clinical_indication_column_code"],
                config["clinical_indication_column_name"],
                config["panel_column"],
                config["test_method_column"]
            ]

    return data
