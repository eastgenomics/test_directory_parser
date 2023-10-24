import argparse

from test_directory_parser import utils
from test_directory_parser import test_directory


def main(args):
    cmd = args.cmd

    if cmd == "rare_disease":
        if args.output:
            output = args.output
        else:
            date = utils.get_date()
            output = f"{date}_RD_TD.json"

        hgnc_data = utils.parse_tsv(args.hgnc)
        rd_test_directory = test_directory.TestDirectory(
            args.test_directory, args.config, "rare_disease", hgnc_data
        )
        rd_test_directory.setup_clinical_indications()

        if args.lab_excel:
            lab_df = utils.parse_lab_excel(args.lab_excel)
            rd_test_directory.filter_clinical_indications(lab_df)

        rd_test_directory.output_json(output)

    elif cmd == "cancer":
        print("Parsing of the cancer test directory is not implemented yet")

    else:
        raise Exception(f"'{cmd}' is not a valid option")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    rare_disease_parser = subparsers.add_parser("rare_disease")
    rare_disease_parser.add_argument(
        "test_directory", help="Path to test directory"
    )
    rare_disease_parser.add_argument(
        "-lab_excel", "--lab_excel",
        help=(
            "Excel file containing the clinical indications that the CUH lab "
            "actually handles"
        )
    )
    rare_disease_parser.set_defaults(which="rare_disease")

    cancer_parser = subparsers.add_parser("cancer")
    cancer_parser.add_argument("test_directory", help="Path to test directory")
    cancer_parser.set_defaults(which="cancer")

    parser.add_argument(
        "-c", "--config",
        help="Config file to know which sheet to gather for example"
    )
    parser.add_argument("-hgnc", "--hgnc", help="Path to the hgnc dump")
    parser.add_argument("-o", "--output", help="Output path and name")

    args = parser.parse_args()
    main(args)
