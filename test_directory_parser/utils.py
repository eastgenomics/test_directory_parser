import datetime

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.schema import MetaData
import numpy as np
import pandas as pd
import regex


def connect_to_local_database(user, passwd, db):
    """ Return cursor of database
    Args:
        user (str): Username for the database
        passwd (str): Password for the user
        db (str): Database name to connect to
    Returns:
        tuple: SQLAlchemy session obj, SQLAlchemy meta obj
    """

    try:
        db = create_engine(
            f"mysql://{user}:{passwd}@127.0.0.1/{db}"
        )
    except Exception as e:
        raise e
    else:
        meta = MetaData()
        meta.reflect(bind=db)
        Session = sessionmaker(bind=db)
        session = Session()
        return session, meta


def get_date():
    """ Return date as string in the following format: YYMMDD

    Returns:
        str: Datetime in string format
    """

    return datetime.datetime.now().strftime("%y%m%d")


def parse_tsv(tsv):
    df = pd.read_csv(tsv, delimiter="\t")
    # replace NA by None so that we can handle them
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
        return data.iloc[0]


def handle_list_panels(panels, hgnc_dump, r_code):
    """ Given a list of "panels", get the hgnc ids/rescue comma panelapp panels

    Args:
        panels (list): List of panels
        hgnc_dump (pd.Dataframe): Hgnc dump dataframe
        r_code (str): R code

    Returns:
        list: List of hgnc ids/panel
    """

    hgnc_ids = []

    # check if the list only contains genes
    # check that the first element has a gene structure i.e. it's
    # not something weird
    if regex.match(r"[A-Z]+[A-Z0-9]+", panels[0]):
        # first element is a gene, so assume that we're dealing with genes
        for panel in panels:
            # for every element in the list, double check that it is a gene
            # because of test directory unpredictableness
            if regex.match(r"[A-Z]+[A-Z0-9]+", panel):
                hgnc_id = find_hgnc_id(panel, hgnc_dump)

                if hgnc_id:
                    hgnc_ids.append(hgnc_id)
                else:
                    # there are instances where the targets can include "and"
                    # at the end of a gene list.
                    # try and catch those instances
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
                        # we didn't manage to find a HGNC id for the gene
                        # symbol
                        print(
                            f"Couldn't find a HGNC id for '{panel}', please "
                            "check manually"
                        )
                        hgnc_ids.append(None)
            else:
                # that element of the list is not a gene
                print(
                    f"This element '{panel}' was not detected as being a gene."
                    " Please check manually"
                )
                hgnc_ids.append(None)

        return hgnc_ids
    else:
        # do regex to see if every element in the list is a panel
        matches = [
            regex.match(r"[A-Za-z0-9-()\ ,]*\([0-9&\ ]+\)", panel)
            for panel in panels
        ]

        if all(matches):
            # all the elements are panelapp panels
            return extract_panelapp_id(panels)
        else:
            # working on it, i realised that trying to rescue the potential
            # panelapp panels is not trivial at all i.e. what if the panel has
            # multiple commas
            # so i think it's safer to manual rescue them
            print(
                f"'{r_code}' has the following potential panelapp panels "
                f"with commas '{panels}'"
            )
            return None


def extract_panelapp_id(panels):
    """ Extract the panelapp id from the target in the test directory

    Args:
        panels (iter): Iterable containing the panels to look at

    Returns:
        list: List of panelapp ids
    """

    panelapp_ids = []

    for panel in panels:
        # find the panelapp id in the panel name
        match = regex.findall(r"\([0-9]+\)", panel)

        if match:
            for m in match:
                # get the whole result and strip out the parentheses
                result = m.strip("()")
                panelapp_ids.append(result)
        else:
            # it might be something like "(panelapp id & panelapp id)", true
            # story
            match = regex.search(r"\([0-9& ]+\)", panel)

            if match:
                result = match.group(0).strip("()")
                panelapp_ids.extend(x.strip() for x in result.split("&"))
            else:
                raise Exception(
                    f"'{panels}' hasn't been recognised as a panelapp panel"
                )

    return panelapp_ids
