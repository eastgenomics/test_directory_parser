import datetime

import numpy as np
import pandas as pd
import regex
from termcolor import cprint


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


def find_hgnc_id(gene_symbol, hgnc_dump):
    # pattern is the gene symbol and only the gene symbol
    pattern = fr"^{gene_symbol}$"
    # try and match the gene symbol in the Approved symbol column
    data = hgnc_dump[
        hgnc_dump["Approved symbol"].str.match(pattern, na=False)
    ]["HGNC ID"]
    # prepare dataframes for aliases and previous symbols
    previous_symbols = hgnc_dump[["Previous symbols", "HGNC ID"]]
    alias_symbols = hgnc_dump[["Alias symbols", "HGNC ID"]]

    # match failed, need to use previous and alias symbols
    if len(data.index) == 0:
        if match_target(r"[A-Z]+[A-Z0-9]+", gene_symbol):
            # previous and alias symbols in hgnc are separated by commas which
            # breaks the pattern matching so i split the appropriate columns
            # creating a df with one column per splitted element
            splitted_previous_symbols = previous_symbols["Previous symbols"].str.split(",", expand=True)
            splitted_alias_symbols = alias_symbols["Alias symbols"].str.split(",", expand=True)
            # copy pasted thing to do the pattern matching for every element of
            # the splitted column 
            previous_symbols_match = np.column_stack(
                [
                    splitted_previous_symbols[col].str.strip().str.match(pattern, na=False)
                    for col in splitted_previous_symbols
                ]
            )
            alias_symbols_match = np.column_stack(
                [
                    splitted_alias_symbols[col].str.strip().str.match(pattern, na=False)
                    for col in splitted_alias_symbols
                ]
            )
            # go back to a dataframe using the numpy array to get the matching
            # rows
            df_previous_symbols = previous_symbols.loc[previous_symbols_match.any(axis=1)]
            df_alias_symbols = alias_symbols.loc[alias_symbols_match.any(axis=1)]

            # check the size of the dataframes and do the appropriate action
            if len(df_previous_symbols) == 0 and len(df_alias_symbols) == 0:
                return None
            elif len(df_previous_symbols) == 1 and len(df_alias_symbols) == 0:
                return df_previous_symbols["HGNC ID"].to_list()[0]
            elif len(df_previous_symbols) == 0 and len(df_alias_symbols) == 1:
                return df_alias_symbols["HGNC ID"].to_list()[0]
            elif len(df_previous_symbols) >= 1 and len(df_alias_symbols) >= 1:
                pass

        # some panel escaped checks
        else:
            raise Exception(f"'{gene_symbol}' escaped checks")
    else:
        return [hgnc_id for hgnc_id in data][0]


def match_target(pattern, targets):
    match = regex.match(pattern, targets)

    if match:
        return True
    else:
        return False


def handle_list_panels(panels, hgnc_dump):
    hgnc_ids = []

    # check if the list only contains genes
    # check that the first element has a gene structure i.e. it's
    # not something weird
    if match_target(r"[A-Z]+[A-Z0-9]+", panels[0]):
        # we are dealing with genes
        for panel in panels:
            if match_target(r"[A-Z]+[A-Z0-9]+", panel):
                hgnc_id = find_hgnc_id(panel, hgnc_dump)

                if hgnc_id:
                    hgnc_ids.append(hgnc_id)
                else:
                    attempt_rescue_gene = [
                        gene.strip() for gene in panel.split("and")
                    ]

                    # check if there were 2 genes from each side of the 'and'
                    if len(attempt_rescue_gene) >= 2:
                        for gene in attempt_rescue_gene:
                            rescued_gene = find_hgnc_id(gene, hgnc_dump)

                            if rescued_gene:
                                hgnc_ids.append(rescued_gene)

        return hgnc_ids
    else:
        # some panelapp panels have commas
        if match_target(r"[A-Za-z0-9-()\ ,]*\([0-9&\ ]+\)", ", ".join(panels)):
            return [", ".join(panels)]

        return None
