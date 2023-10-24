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


def parse_lab_excel(excel_file):
    df = pd.read_excel(excel_file)
    return df[df["NGS Technology"].isin(["CEN", "WES"])]


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

    df_res = pd.Series(
        [gene_symbol, None, None, None], index=[
            "Gene symbol", "HGNC ID", "Previous", "Alias"
        ]
    )

    hgnc_id = None

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
                return df_res

            elif len(df_previous_symbols) == 1 and len(df_alias_symbols) == 0:
                # found only a previous symbol, return the HGNC id
                hgnc_id = df_previous_symbols["HGNC ID"].to_list()[0]
                df_res.at["Previous"] = True
                df_res.at["Alias"] = False

            elif len(df_previous_symbols) == 0 and len(df_alias_symbols) == 1:
                # found only a alias symbol, return the HGNC id
                hgnc_id = df_alias_symbols["HGNC ID"].to_list()[0]
                df_res.at["Previous"] = False
                df_res.at["Alias"] = True

            elif len(df_previous_symbols) >= 1 and len(df_alias_symbols) >= 1:
                # found previous and alias symbols, cry
                df_res.at["Previous"] = True
                df_res.at["Alias"] = True

            elif len(df_previous_symbols) >= 1:
                df_res.at["Previous"] = True
                df_res.at["Alias"] = False

            elif len(df_alias_symbols) >= 1:
                df_res.at["Previous"] = False
                df_res.at["Alias"] = True

    else:
        hgnc_id = data.iloc[0]

    df_res.at["HGNC ID"] = hgnc_id

    return df_res
