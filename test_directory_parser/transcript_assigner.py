import datetime

from test_directory_parser import utils


def get_date_for_db():
    """ Return date in YYYY-MM-DD format for insertion in database

    Returns:
        str: String with date
    """

    return datetime.datetime.now().strftime("%Y-%m-%d")


def parse_mane_file(mane_file, hgnc_dump):
    """ Parse a MANE CSV file downloaded from
    http://tark.ensembl.org/web/mane_GRCh37_list/

    Args:
        mane_file (str): Path to MANE CSV file
        hgnc_dump (dict): Dict containing the HGNC data

    Returns:
        dict: Dict containing MANE data
    """

    mane_data = {}

    with open(mane_file) as f:
        for line in f:
            # each column has quotes around them so i'm striping them i.e.
            # "A1BG","MANE SELECT","ENST00000263100.8","NM_130786.4"
            cleaned_data = [ele.strip("\"") for ele in line.strip().split(",")]
            gene_symbol, mane_type, ensembl_38, refseq = cleaned_data[0:4]

            if mane_type == "MANE SELECT":
                hgnc_id = utils.find_hgnc_id(gene_symbol, hgnc_dump)

                if hgnc_id:
                    mane_data[hgnc_id] = refseq
                else:
                    print(f"Couldn't find hgnc id for {gene_symbol}")

    return mane_data


def parse_gene_transcript(g2t_file):
    """ Parse gene to transcripts file (2 columns file) for missing genes in
    the panel database

    Args:
        g2t_file (str): Path to the g2t file

    Returns:
        dict: Dict with gene as the key and the transcripts as values
    """

    data = {}

    with open(g2t_file) as f:
        for line in f:
            hgnc_id, transcript = line.strip().split()
            data.setdefault(hgnc_id, []).append(transcript)

    return data


def find_HGMD_transcript(session, meta, hgnc_id):
    """ Find the HGMD transcript using a HGNC id

    Args:
        session (session object): Session SQLAlchemy object
        meta (meta object): Meta SQLAlchemy object
        hgnc_id (str): HGNC id

    Raises:
        Exception: if 2 or more transcripts are found for the given HGNC id

    Returns:
        str: HGMD transcript
    """

    markgene_tb = meta.tables["markname"]
    gene2refseq_tb = meta.tables["gene2refseq"]

    # get the hgmd transcripts using the HGNC id provided
    hgmd_transcripts = session.query(
        gene2refseq_tb.c.refcore, gene2refseq_tb.c.refversion
    ).join(
        markgene_tb, markgene_tb.c.gene_id == gene2refseq_tb.c.hgmdID
    ).filter(markgene_tb.c.hgncID == hgnc_id[5:]).all()

    hgmd_transcript = [
        f"{hgmd_base}.{hgmd_version}"
        for hgmd_base, hgmd_version in hgmd_transcripts
    ]

    if len(hgmd_transcript) == 0:
        return
    elif len(hgmd_transcript) == 1:
        return hgmd_transcript[0]
    else:
        # HGNC ids should only be linked to one HGMD transcript
        raise Exception(
            f"{hgnc_id} has the following transcripts in the HGMD database: "
            f"{', '.join(hgmd_transcript)}"
        )


def assign_transcripts(session, meta, mane_select_data, g2t_data):
    """ Assign a clinical transcript status for all the genes in the g2t file

    Args:
        session (session object): Session SQLAlchemy object
        meta (meta object): Meta SQLAlchemy object
        mane_select_data (dict): Dict containing the MANE select data
        g2t_data (dict): Dict containing the g2t data

    Returns:
        dict: Dict of dict with the gene as key, status as subkey and the
        transcripts as values
    """

    data = {}

    for gene, transcripts in g2t_data.items():
        data.setdefault(gene, {})
        data[gene].setdefault("no_clinical_transcript", [])

        for tx in transcripts:
            tx_base, tx_version = tx.split(".")

            if gene in mane_select_data:
                mane_transcript = mane_select_data[gene]

                mane_base, mane_version = mane_transcript.split(".")

                # compare transcripts without the versions
                if tx_base == mane_base:
                    data[gene]["clinical_transcript"] = [tx, "MANE"]
                    continue

            # if we already have a clinical transcript, that means that we
            # already have a MANE transcript
            if "clinical_transcript" in data[gene]:
                data[gene]["no_clinical_transcript"].append(tx)
                continue

            hgmd_transcript = find_HGMD_transcript(session, meta, gene)

            if hgmd_transcript:
                hgmd_base, hgmd_version = hgmd_transcript.split(".")

                if tx_base == hgmd_base:
                    data[gene]["clinical_transcript"] = [tx, "HGMD"]
                    continue

            data[gene]["no_clinical_transcript"].append(tx)

    return data


def write_sql_queries(data, output_path):
    """ Write SQL queries to run on the panel database to add genes and
    transcripts output

    Args:
        data (dict): Dict of dict with the gene as key, status as subkey and
        the transcripts as values 
        output_path (str): Path to output file
    """

    with open(output_path, "w") as f:
        for gene in data:
            f.write(f"INSERT INTO gene (hgnc_id) VALUES (\"{gene}\");\n")
            f.write("SET @gene_id = (SELECT LAST_INSERT_ID());\n")

            for status, txs in data[gene].items():
                if status == "clinical_transcript":
                    tx, source = txs
                    tx_base, tx_version = tx.split(".")
                    f.write((
                        "INSERT INTO transcript (transcript, version, canonical) "
                        f"VALUES (\"{tx_base}\", \"{tx_version}\", 0);\n"
                    ))
                    f.write("SET @transcript_id = (SELECT LAST_INSERT_ID());\n")
                    f.write((
                        "INSERT INTO genes2transcripts "
                        "(clinical_transcript, date, gene_id, reference_id, transcript_id) "
                        f"VALUES (1, \"{get_date_for_db()}\", @gene_id, 1, @transcript_id);\n"
                    ))

                else:
                    for tx in txs:
                        tx_base, tx_version = tx.split(".")
                        f.write((
                            "INSERT INTO transcript (transcript, version, canonical) "
                            f"VALUES (\"{tx_base}\", \"{tx_version}\", 0);\n"
                        ))
                        f.write("SET @transcript_id = (SELECT LAST_INSERT_ID());\n")
                        f.write((
                            "INSERT INTO genes2transcripts "
                            "(clinical_transcript, date, gene_id, reference_id, transcript_id) "
                            f"VALUES (0, \"{get_date_for_db()}\", @gene_id, 1, @transcript_id);\n"
                        ))


def write_transcript_status(data):
    """ Write 2 files with genes depending on their clinical status

    Args:
        data (dict): Dict of dict with the gene as key, status as subkey and
        the transcripts as values 
    """

    f1 = open("genes_with_clinical_transcripts.txt", "w")
    f2 = open("genes_with_no_clinical_transcripts.txt", "w")

    for gene, transcript_data in data.items():
        # check if a gene has "clinical_transcript" subkey
        all_clinical_transcripts = any([
            True if status == "clinical_transcript" else False
            for status in transcript_data
        ])

        if all_clinical_transcripts:
            f1.write(f"{gene}\n")
        else:
            f2.write(f"{gene}\n")

    f1.close()
    f2.close()
