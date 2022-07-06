import datetime

import numpy as np
import pandas as pd
import regex


def get_date():
    """ Return date as string in the following format: YYMMDD

    Returns:
        str: Datetime in string format
    """

    return datetime.datetime.now().strftime("%y%m%d")


def parse_tsv(tsv):
    df = pd.read_csv(tsv, delimiter="\t")
    df_with_none = df.where(pd.notnull(df), None)
    return df_with_none


def find_hgnc_id(gene_symbol, hgnc_dump):
    """ Find hgnc id using the hgnc dump

    Args:
        gene_symbol (str): Gene symbol
        hgnc_dump (pd.Dataframe): Hgnc dump dataframe

    Raises:
        Exception: if a panel has escaped previous checks

    Returns:
        str: Hgnc id
    """

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
        if regex.match(r"[A-Z]+[A-Z0-9]+", gene_symbol):
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
                # couldn't find a previous or alias symbol
                return None
            elif len(df_previous_symbols) == 1 and len(df_alias_symbols) == 0:
                # found only a previous symbol, return the HGNC id
                return df_previous_symbols["HGNC ID"].to_list()[0]
            elif len(df_previous_symbols) == 0 and len(df_alias_symbols) == 1:
                # found only a alias symbol, return the HGNC id
                return df_alias_symbols["HGNC ID"].to_list()[0]
            elif len(df_previous_symbols) >= 1 and len(df_alias_symbols) >= 1:
                # found previous and alias symbols, cry
                print(
                    "Couldn't find a non ambiguous HGNC id for "
                    f"'{gene_symbol}'"
                )

        # some panel escaped checks
        else:
            raise Exception(f"'{gene_symbol}' escaped checks")
    else:
        return [hgnc_id for hgnc_id in data][0]


def handle_list_panels(panels, hgnc_dump):
    """ Given a list of "panels", get the hgnc ids/rescue comma panelapp panels

    Args:
        panels (list): List of panels
        hgnc_dump (pd.Dataframe): Hgnc dump dataframe

    Returns:
        list: List of hgnc ids/panel
    """

    hgnc_ids = []

    # check if the list only contains genes
    # check that the first element has a gene structure i.e. it's
    # not something weird
    if regex.match(r"[A-Z]+[A-Z0-9]+", panels[0]):
        # we are dealing with genes
        for panel in panels:
            if regex.match(r"[A-Z]+[A-Z0-9]+", panel):
                hgnc_id = find_hgnc_id(panel, hgnc_dump)

                if hgnc_id:
                    hgnc_ids.append(hgnc_id)
                else:
                    if "and" in panel:
                        attempt_rescue_gene = [
                            gene.strip() for gene in panel.split("and")
                        ]

                        # check if there were 2 genes from each side of the 'and'
                        if len(attempt_rescue_gene) >= 2:
                            for gene in attempt_rescue_gene:
                                rescued_gene = find_hgnc_id(gene, hgnc_dump)

                                if rescued_gene:
                                    hgnc_ids.append(rescued_gene)
                                else:
                                    hgnc_ids.append(None)
                    else:
                        hgnc_ids.append(None)

        return hgnc_ids
    else:
        # some panelapp panels have commas
        if regex.match(r"[A-Za-z0-9-()\ ,]*\([0-9&\ ]+\)", ", ".join(panels)):
            return [", ".join(panels)]

        return None
