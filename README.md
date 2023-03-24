# test_directory_parser
Script(s) to parse the test directories

## Before running the code

### HGNC dump

To generate the HGNC dump, you can go to: https://www.genenames.org/download/custom/

And for the code to work, the following checkboxes need to be checked when you download the dump:

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
    "name": "220421_RD",
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

Setup your environment first:

```bash
python3 -m venv ${path_to_env}/${env_name}
source ${path_to_env}/${env_name}/bin/activate
pip install -r requirements.txt
```

There are 2 main modes for this script:

- Parse the test directory (only rare disease for now)

```bash
# outputs a json containing cleaned data from the given test directory
python main.py --hgnc ${hgnc_dump.txt} rare_disease ${test_directory.xslx} configs/${config} [-o ${output_path}]
```

- Check content of given test directory against panel database and output:
  - absent_genes_per_panel.txt: TSV with clinical indications and the potential genes that are absent from the database
  - no_clinical_transcript_per_clinical_indication.txt: TSV with clinical indications and the genes that don't have clinical transcripts
  - absent_genes_per_panel.txt: Text file with all genes that need to be added in the database (removes some duplicates)
  - genes_with_no_clinical_transcripts.txt: Text file with all genes that have no clinical transcripts in the database
  - (optional) clinical_indication_with_${string_filter}.txt: Text file with clinical indication codes filtered on the change column using a string filter

```bash
# output mandatory check files
python main.py checker ${panel_database_username} ${panel_database_passwd} configs/${config}
# output all the files
python main.py checker ${panel_database_username} ${panel_database_passwd} configs/${config} -f ${filter_string}
```
