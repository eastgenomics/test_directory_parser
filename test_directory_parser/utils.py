import datetime

import pandas as pd
import regex


def get_date():
    """ Return date as string in the following format: YYMMDD

    Returns:
        str: Datetime in string format
    """

    now = datetime.datetime.now()
    date = now.strftime("%y%m%d")
    return date


def parse_tsv(tsv):
    data = pd.read_csv(tsv, delimiter="\t")
    noned_data = data.where(pd.notnull(data), None)
    return noned_data


def match_target(pattern, targets):
    match = regex.match(pattern, targets)

    if match:
        return True
    else:
        return False
