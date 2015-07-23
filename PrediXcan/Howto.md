# INTRODUCTION

*predict_gene_expression.py* is a command line tool to,
 as you may imagine, generate the gene expression matrix.
 
You need Python and numpy installed as prerequisites.

# Input Data

This script needs the following data:

- genotype biallelic SNP dosage data
- a Transcriptome Prediction Model (see sample files below)

It expects this input in a folder structure (at the same level the script is run) such as:

```bash
data/dosages     # folder containing gzipped dosages, such as "chr1.dosage.txt.gz"
data/weights.db  # sqlite database with weights from transcriptome model
```

If the input files are named exactly as in the previous example, you can just run:

```bash
$ ./predict_gene_expression.py
```

And it will output a txt file with the transcription matrix.

## Parameters

The script supports some command line arguments. Type:
```bash
$ ./predict_gene_expression.py --help
```
for reference. 

For example:

```bash
$ ./predict_gene_expression.py \
    --weights data/DGN-WB_0.5.db \
    --dosages another_folder \
    --output results/expression.txt
```

will cause the dosages at *another_folder* and the database *data/DGN-WB_0.5.db* to be read,
and the result will be saved in *results/expression.txt*

## Sample data

### Transcriptome models

You can download sqlite database models from [here](https://app.box.com/s/5nejbvzgsis77wtrxt8xokyn7pcdgnnp)

### Dosage data

This files are expected to be compressed files in the following format:
```
# chr1.dosage.txt.gz sample content, made up of the following columns:
# CHROMOSOME RSID POSITION ALLELE_1 ALLELE_2 ALLELE_2_FREQUENCY_AVERAGE (...)
# where (...) means a list of allele frequency values for different persons
...
chr1 rs10399749 55299 C T 0.164302600472813 0 0 1 ...
...
```

If you have access to Hae Kyung Im's public data repository, and have Amazon Web Services command line tools,
you can get dosage by executing:

``` bash
mkdir data
mkdir data/dosages
cd data/dosages
aws s3 cp s3://imlab-open/Data/1000Genomes/Transcriptome/GEUVADIS/dosagefiles-hapmap2/ . --recursive
```