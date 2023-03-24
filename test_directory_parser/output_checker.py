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

    latest_version = max(
        [Version(version.split("|")[0]) for version in versions]
    )

    multiple_latest_versions = [
        version for version in versions if str(latest_version) in version
    ]

    if len(multiple_latest_versions) == 1:
        return multiple_latest_versions[0]

    elif len(multiple_latest_versions) == 0:
        raise f"Couldn't get the latest version given: '{versions}'"

    else:
        add_on_versions = [
            version for version in multiple_latest_versions if "|" in version
        ]

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
            raise (
                "No add on versions detected where there shouldn't be: "
                f"'{versions}'"
            )


def get_current_panel_genes(session, meta, panelapp_id):
    panel_tb = meta.tables["panel"]
    feature_tb = meta.tables["feature"]
    panel_features_tb = meta.tables["panel_features"]
    gene_tb = meta.tables["gene"]

    res = session.query(gene_tb.c.hgnc_id, panel_features_tb.c.panel_version)\
        .select_from(gene_tb)\
        .join(feature_tb).join(panel_features_tb).join(panel_tb)\
        .filter(panel_tb.c.panelapp_id == panelapp_id).distinct().all()

    genes = {}

    for hgnc_id, panel_version in res:
        genes.setdefault(panel_version, []).append(hgnc_id)

    if len(genes) == 1:
        return set(list(genes.values())[0])

    elif len(genes) >= 2:
        latest_version = extract_latest_version(genes.keys())
        return genes[str(latest_version)]

    else:
        return


def check_genes_in_database(session, meta, genes):
    gene_tb = meta.tables["gene"]
    g2t_tb = meta.tables["genes2transcripts"]

    absent_genes = set()
    no_clinical_transcripts = set()

    for gene in genes:
        query = session.query(gene_tb.c.hgnc_id, g2t_tb.c.clinical_transcript)\
            .join(g2t_tb)\
            .filter(gene_tb.c.hgnc_id == gene).all()

        if query:
            has_clinical_tx = False
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

    panel_absent_genes = {}
    panel_no_clinical_transcripts = {}

    all_absent_genes = set()
    all_no_clinical_transcripts = set()

    for indication in td_parser_output["indications"]:
        panels = indication["panels"]
        code = indication["code"]

        if panels:
            for panelapp_id in panels:
                if regex.search(r"^[0-9]+", panelapp_id):
                    current_genes = get_current_panel_genes(
                        session, meta, panelapp_id
                    )
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


def write_dict(data_dict, name):
    with open(name, "w") as f:
        for key, values in data_dict.items():
            for value in values:
                f.write(f"{key}\t{value}\n")


def write_list(iterable_of_things, name):
    with open(name, "w") as f:
        for gene in iterable_of_things:
            f.write(f"{gene}\n")


def get_clinical_indications(td_parser_output, filter_option):
    output = []

    for indication in td_parser_output["indications"]:
        if regex.search(filter_option, indication["changes"], re.IGNORECASE):
            output.append(indication["code"])

    return output
