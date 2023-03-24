import re

from packaging.version import Version
import regex


def extract_latest_version(versions):
    """ Given an iterable of versions, returns the latest version

    Args:
        versions (iter): Iterable of version parsable strings

    Returns:
        Version: Version type object
    """

    # get the latest version using the base version since we might have a
    # mixture of add on panels (2.0|1) and normal panels (2.0, 2.1)
    # i.e. ["2.0", "2.1", "2.1|1"] --> "2.1"
    latest_version = max(
        [Version(version.split("|")[0]) for version in versions]
    )

    # since we have add on version we need to get all versions that contain
    # the previously gathered max version
    # i.e. ["2.0", "2.1", "2.1|1"] --> ["2.1", "2.1|1"]
    multiple_latest_versions = [
        version for version in versions if str(latest_version) in version
    ]

    # if we get only 1 element in the list, we can confidently return the max
    # version
    if len(multiple_latest_versions) == 1:
        return multiple_latest_versions[0]

    elif len(multiple_latest_versions) == 0:
        raise f"Couldn't get the latest version given: '{versions}'"

    else:
        # if we have multiple elements in the list, that probably means that
        # we have at least one add on version in the list.
        # We need to check whether we have only one or multiple
        add_on_versions = [
            version for version in multiple_latest_versions if "|" in version
        ]

        # if we have multiple add on, do a version check on the version after
        # the pipe to finally return the latest of the add on panel versions
        if len(add_on_versions) >= 1:
            max_add_on_version = max([
                Version(version.split("|")[1])
                for version in add_on_versions
            ])

            return [
                version for version in versions
                if f"|{str(max_add_on_version)}" in version
            ][0]

        else:
            # i mean, i wouldn't understand what is going on there
            raise (
                "No add on versions detected where there shouldn't be: "
                f"'{versions}'"
            )


def get_current_panel_genes(session, meta, panelapp_id):
    """From a panelapp id, get the content of the panel from the panel database

    Args:
        session (Session object): SQLAlchemy Session object connected to the
            panel database
        meta (Meta object): SQLAlchemy Meta object
        panelapp_id (str): String with the panelapp id of interest

    Returns:
        set: Set containing genes of the latest panel version of the panel
    """

    panel_tb = meta.tables["panel"]
    feature_tb = meta.tables["feature"]
    panel_features_tb = meta.tables["panel_features"]
    gene_tb = meta.tables["gene"]

    # query the database for the hgnc ids and panel version of a panel
    # using its panelapp id
    res = session.query(gene_tb.c.hgnc_id, panel_features_tb.c.panel_version)\
        .select_from(gene_tb)\
        .join(feature_tb).join(panel_features_tb).join(panel_tb)\
        .filter(panel_tb.c.panelapp_id == panelapp_id).distinct().all()

    genes = {}

    # for each row of result, store the genes using the panel version as a key
    for hgnc_id, panel_version in res:
        genes.setdefault(panel_version, []).append(hgnc_id)

    # only 1 panel version recovered, return the values
    if len(genes) == 1:
        return set(list(genes.values())[0])

    # 2 or more version detected, so get latest version
    elif len(genes) >= 2:
        latest_version = extract_latest_version(genes.keys())
        return genes[str(latest_version)]

    else:
        return


def check_genes_in_database(session, meta, genes):
    """ Check if the genes given are present in the database and if they have
    clinical transcripts

    Args:
        session (Session object): SQLAlchemy Session object connected to the
            panel database
        meta (Meta object): SQLAlchemy Meta object
        genes (iterable): Iterable of genes (HGNC ids)

    Returns:
        list: List of sets for genes not present in the database and genes
        that don't have clinical transcripts
    """

    gene_tb = meta.tables["gene"]
    g2t_tb = meta.tables["genes2transcripts"]

    absent_genes = set()
    no_clinical_transcripts = set()

    for gene in genes:
        # query the database using the HGNC id and join the gene and g2t tables
        query = session.query(gene_tb.c.hgnc_id, g2t_tb.c.clinical_transcript)\
            .join(g2t_tb)\
            .filter(gene_tb.c.hgnc_id == gene).all()

        if query:
            has_clinical_tx = False

            # for every hgnc id/transcript couple, check if the clinical
            # transcript status == 1
            # i.e. the transcript is the clinical transcript
            for hgnc_id, clinical_transcript in query:
                if clinical_transcript == 1:
                    has_clinical_tx = True

            if not has_clinical_tx:
                no_clinical_transcripts.add(gene)

        else:
            absent_genes.add(gene)

    return absent_genes, no_clinical_transcripts


def compare_panelapp_panels_content(
    session, meta, td_parser_output, signedoff_panels
):
    """ Compare the content of the latest signedoff version of a panelapp panel
    with what is present in the database

    Args:
        session (Session object): SQLAlchemy Session object connected to the
            panel database
        meta (Meta object): SQLAlchemy Meta object
        td_parser_output (dict): Dict of the json td parser output
        signedoff_panels (dict): Dict of the signedoff panelapp panels
        (key = panelapp id, value = Panelapp object)

    Returns:
        list: List of iterables:
            - dict: absent genes per panel
            - dict: genes with no clinical transcript per panel
            - set: absent genes total
            - set: genes with no clinical transcript total
    """

    panel_absent_genes = {}
    panel_no_clinical_transcripts = {}

    all_absent_genes = set()
    all_no_clinical_transcripts = set()

    # loop the clinical indication of the test directory
    for indication in td_parser_output["indications"]:
        panels = indication["panels"]
        code = indication["code"]

        if panels:
            # clinical indications can have multiple panels
            for panelapp_id in panels:
                # check if the target is a panelapp panel or single genes
                if regex.search(r"^[0-9]+", panelapp_id):
                    # get the genes in the database
                    current_genes = get_current_panel_genes(
                        session, meta, panelapp_id
                    )
                    # get the green genes of panelapp panel
                    panelapp_genes = signedoff_panels[int(panelapp_id)]\
                        .get_genes()
                    panelapp_genes = set([
                        gene["hgnc_id"] for gene in panelapp_genes
                    ])

                    if not current_genes:
                        print(
                            f"{code} might be missing from the current database"
                        )
                        continue

                    # this will give us the new genes coming in the update
                    new_genes = panelapp_genes.difference(current_genes)

                    (
                        absent_genes, no_clinical_transcripts
                    ) = check_genes_in_database(session, meta, new_genes)

                    panel_absent_genes.setdefault(code, set())\
                        .update(absent_genes)
                    panel_no_clinical_transcripts.setdefault(code, set())\
                        .update(no_clinical_transcripts)

                    all_absent_genes.update(absent_genes)
                    all_no_clinical_transcripts.update(
                        no_clinical_transcripts
                    )

    return (
        panel_absent_genes, panel_no_clinical_transcripts,
        all_absent_genes, all_no_clinical_transcripts
    )


def write_data(data, name):
    with open(name, "w") as f:
        if isinstance(data, dict):
            for key, values in data.items():
                for value in values:
                    f.write(f"{key}\t{value}\n")

        elif isinstance(data, (list, tuple, set)):
            for d in data:
                f.write(f"{d}\n")


def get_clinical_indications(td_parser_output, filter_option):
    output = []

    for indication in td_parser_output["indications"]:
        if regex.search(filter_option, indication["changes"], re.IGNORECASE):
            output.append(indication["code"])

    return output
