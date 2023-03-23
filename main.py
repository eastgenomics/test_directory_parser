import argparse
import json

from panelapp import queries

from test_directory_parser import rare_disease
from test_directory_parser import utils
from test_directory_parser import test_directory
from test_directory_parser import output_checker


def main(args):
    cmd = args.cmd

    if cmd == "rare_disease":
        if args.output:
            output = args.output
        else:
            date = utils.get_date()
            output = f"{date}_RD_TD.json"

        hgnc_data = utils.parse_tsv(args.hgnc)
        config_data = rare_disease.parse_config(args.config)
        sheet, change_column = rare_disease.parse_rare_disease_td(
            args.test_directory, config_data
        )
        rd_test_directory = test_directory.TestDirectory(
            sheet, change_column, config_data, "rare_disease", hgnc_data
        )
        rd_test_directory.setup_clinical_indications()
        rd_test_directory.output_json(output)

    elif cmd == "cancer":
        print("Parsing of the cancer test directory is not implemented yet")

    elif cmd == "checker":
        session, meta = utils.connect_to_panel_database(
            args.username, args.passwd
        )
        td_parser_output = json.load(open(args.td_parser_output))

        cis = output_checker.get_clinical_indications(
            td_parser_output, args.filter
        )

        panelapp_panels = queries.get_all_signedoff_panels()

        (
            absent_genes_per_panel, no_clinical_transcripts_per_panel,
            absent_genes, genes_no_clinical_transcripts
        ) = output_checker.compare_panelapp_panels_content(
            session, meta, td_parser_output, panelapp_panels
        )

        output_checker.write_dict(
            absent_genes_per_panel, "absent_genes_per_panel.txt"
        )
        output_checker.write_dict(
            no_clinical_transcripts_per_panel,
            "no_clinical_transcritps_per_panel.txt"
        )
        output_checker.write_list(
            cis, f"clinical_indication_with_{args.filter}.txt"
        )
        output_checker.write_list(absent_genes, "genes_not_in_database.txt")
        output_checker.write_list(
            genes_no_clinical_transcripts,
            "genes_with_no_clinical_transcripts.txt"
        )
    else:
        raise Exception(f"'{cmd}' is not a valid option")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    rare_disease_parser = subparsers.add_parser("rare_disease")
    rare_disease_parser.add_argument(
        "test_directory", help="Path to test directory"
    )
    rare_disease_parser.set_defaults(which="rare_disease")

    cancer_parser = subparsers.add_parser("cancer")
    cancer_parser.add_argument("test_directory", help="Path to test directory")
    cancer_parser.set_defaults(which="cancer")

    checker_parser = subparsers.add_parser("checker")
    checker_parser.add_argument(
        "username", help="Username for the panel database"
    )
    checker_parser.add_argument(
        "passwd", help="Password for given username"
    )
    checker_parser.add_argument(
        "td_parser_output", help="Json output of the test directory parser"
    )
    checker_parser.add_argument(
        "-f", "--filter", help="Filter option on the change column"
    )
    checker_parser.set_defaults(which="checker")

    parser.add_argument(
        "-c", "--config",
        help="Config file to know which sheet to gather for example"
    )
    parser.add_argument("-hgnc", "--hgnc", help="Path to the hgnc dump")
    parser.add_argument("-o", "--output", help="Output path and name")

    args = parser.parse_args()
    main(args)
