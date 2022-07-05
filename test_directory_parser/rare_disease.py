import json

import pandas as pd


def parse_config(config):
    """ Parse config file

    Args:
        config (str): Path to the config file

    Returns:
        dict: Dict of config data
    """

    config_content = open(config)
    config_data = json.load(config_content)
    config_content.close()
    return config_data


def parse_rare_disease_td(test_directory, config):
    """Parse rare disease test directory using the config file

    Args:
        test_directory (str): Path to the test directory
        config (dict): Dict containing the data for the config file

    Returns:
        pandas.Dataframe: Dataframe containing the columns of interest
    """

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
