import argparse

from test_directory_parser import rare_disease
from test_directory_parser import utils
from test_directory_parser import test_directory


def main(args):
    cmd = args.cmd

    if cmd == "rare_disease":
        if args.output:
            output = args.output
        else:
            date = utils.get_date()
            output = f"{date}_RD_TD_output"

        hgnc_data = utils.parse_tsv(args.hgnc)
        config_data = rare_disease.parse_config(args.config)
        sheet = rare_disease.parse_rare_disease_td(
            args.test_directory, config_data
        )
        rd_test_directory = test_directory.TestDirectory(
            sheet, config_data, "rare_disease", hgnc_data
        )
        rd_test_directory.setup_clinical_indications()
        rd_test_directory.output()

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
    rare_disease_parser.set_defaults(which="rare_disease")

    cancer_parser = subparsers.add_parser("cancer")
    cancer_parser.add_argument("test_directory", help="Path to test directory")
    cancer_parser.set_defaults(which="cancer")

    parser.add_argument(
        "config", help="Config file to know which sheet to gather for example"
    )
    parser.add_argument("-hgnc", "--hgnc", help="Path to the hgnc dump")
    parser.add_argument("-o", "--output", help="Output path and name")

    args = parser.parse_args()
    main(args)
