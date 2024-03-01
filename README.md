# test_directory_parser
Script(s) to parse the test directories.

These test directories have gone through a checking phase using https://github.com/eastgenomics/test_directory_checker

## Before running the code

### HGNC dump

To generate the HGNC dump, you can go to: https://www.genenames.org/download/custom/

And for the code to work, the following checkboxes need to be checked when you download the dump (they are by default):

- HGNC ID
- Approved symbol
- Alias symbols
- Previous symbols

### Test directory

The rare disease and cancer test directories are downloadable here: https://www.england.nhs.uk/publication/national-genomic-test-directories/

Right now this version of the code can only be used to parse the rare disease test directory.

### Config file

The config file is used to indicate the header line, the name of the sheet of interest and the name of the columns that need to be gathered and processed.

Right now the columns containing the test code, clinical indication name, test methods and the targets columns are processed without addition of code.

The `ngs_test_methods` field contains the test methods that we want to keep, so filtering will be applied on the test method using this list

```json
{
    "name": "220421_RD",
    "sheet_of_interest": "R&ID indications",
    "clinical_indication_column_code": "Test ID",
    "clinical_indication_column_name": "Clinical Indication",
    "panel_column": "Target/Genes",
    "test_method_column": "Test Method",
    "specialty_column": "Specialist test group",
    "header_index": 1,
    "specialism_of_interest": ["Core", "Endocrinology", "Neurology"],
    "ngs_test_methods": [
        "Medium panel", "Single gene sequencing <=10 amplicons",
        "Single gene sequencing <10 amplicons",
        "Single gene sequencing >=10 amplicons",
        "Single gene testing (<10 amplicons)", "small panel", "Small panel",
        "WES or Large panel", "WES or Large Panel", "WES or Large penel",
        "WES or Medium panel", "WES or Medium Panel", "WES or Small Panel", "WGS"
    ]
}
```

### Python environment

Setup your environment first:

```bash
python3 -m venv ${path_to_env}/${env_name}
source ${path_to_env}/${env_name}/bin/activate
pip install -r requirements.txt
```

## How to run

```bash
# outputs a json containing cleaned data from the given test directory
python main.py -c configs/${config} [-o ${output_path}] --hgnc ${hgnc_dump.txt} rare_disease ${test_directory.xlsx} 
```

## Run unittests

```bash
python -m unittest test_directory_parser.tests
```

## Output

The code will output a JSON file with the following default name `${YYMMDD}_RD_TD.json` with the following format:

```json
{
  "td_source": "name_of_td_file_used_at_runtime",
  "config_source": "config_file_named_used",
  "date": "date_at_runtime",
  "indications": [
    {
      "name": "CI1",
      "code": "R1.1",
      "gemini_name": "R1.1_CI1_P",
      "test_method": "test_method1",
      "panels": [
        "panelapp_id"
      ],
      "original_targets": "Panel 1 (panelapp_id)",
      "changes": "No changes"
    },
    {
      "name": "CI2",
      "code": "R2.1",
      "gemini_name": "R2.1_CI2_P",
      "test_method": "test_method2",
      "panels": [
        "HGNC_ID"
      ],
      "original_targets": "Gene symbol",
      "changes": "No changes"
    },
  ]
}
```

This output is than used to import this data into the panel database using panel_ops (https://github.com/eastgenomics/panel_ops).
