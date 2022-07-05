# test_directory_parser
Script(s) to parse the test directories

## Before running the code

### HGNC dump

To generate the HGNC dump, you can go to: https://www.genenames.org/download/custom/

And for the code to work, the following checkboxes need to be checked:

- HGNC ID
- Approved symbol
- Alias symbols
- Previous symbols

Example:

```tsv
HGNC ID	Approved symbol	Previous symbols	Alias symbols
HGNC:1	A12M1		
HGNC:10	A2MRAP		
HGNC:100	ASIC1	ACCN2	BNaC2, hBNaC2
HGNC:1000	BCL5		
HGNC:10000	RGS4	SCZD9	
```

### Test directory

The rare disease and cancer test directories are downloadable here: https://www.england.nhs.uk/publication/national-genomic-test-directories/

Right now this version of the code can only be used to parse the rare disease test directory.

### Config file

The config file is used to indicated the header line, the name of the sheet of interest and the name of the columns that need to be gathered and processed.

Right now the columns containing the test code, clinical indication name, test methods and the targets columns are processed without addition of code.

The `ngs_test_methods` field contains the test methods that we want to keep, so filtering will be applied on the test method using this list

```json
{
    "sheet_of_interest": "R&ID indications",
    "clinical_indication_column_code": "Test ID",
    "clinical_indication_column_name": "Clinical Indication",
    "panel_column": "Target/Genes",
    "test_method_column": "Test Method",
    "header_index": 1,
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

## How to run

```bash
source envs/test_directory_parser/bin/activate

python main.py --hgnc ${hgnc_dump.txt} rare_disease ${test_directory.xlsx} configs/220421_RD.json
```
