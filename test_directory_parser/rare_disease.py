import json

import pandas as pd


def parse_config(config: str):
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


def parse_rare_disease_td(test_directory: str, config: dict):
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
            # find name of change column
            change_column = [
                col for col in xls[sheet].columns
                if config['changes_column'].lower() in col.lower()
            ]

            if change_column:
                if len(change_column) == 1:
                    change_column = change_column[0]
                else:
                    raise Exception((
                        "2 or more columns were detected as having "
                        f"'{config['changes_column']}' "
                        "in their name breaking the changes gathering.\n"
                        f"Matched column names: {';'.join(change_column)}"
                    ))
            else:
                raise Exception("Couldn't find the change column.")

            data = xls[sheet].loc(axis=1)[
                config["clinical_indication_column_code"],
                config["clinical_indication_column_name"],
                config["panel_column"],
                config["test_method_column"],
                config["specialty_column"],
                change_column
            ]

    # filter using the specialisms used in the lab
    filtered_data = data.loc[
        data[
            config["specialty_column"]
        ].isin(config["specialism_of_interest"])
    ]

    return filtered_data, change_column
