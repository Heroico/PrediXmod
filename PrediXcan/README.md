PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 

Reference
---------
PrediXcan: Trait Mapping Using Human Transcriptome Regulation
Eric R. Gamazon†, Heather E. Wheeler†, Sahar Mozaffari, Kaanan P. Shah, Keston Aquino-Michaels, GTEx Consortium, Dan L. Nicolae, Nancy J. Cox, and Hae Kyung Im* (Under revision)

Instructions
------------

You will need 

Input: 

- genotype file 
- phenotype file
- transcriptome prediction model (sqlite db to be downloaded from [here](https://app.box.com/s/i0gxg9703jj9lypq3umopponnsdqupze "Titulo")
 	- tissue: Whole Blood (default)
	- source: DGN (default)
	- model: Elastic Net 0.50 (default)


download this python script to compute predicted expression levels, and the cross validated performance measures for each gene. 

The script predict_gene_expression.py predicts gene expression levels using prediction models (stored in sqlite db such as Skin-SunExposed(Lowerleg)_0.5.db) and whole genome variation data.

For now the association with phenotype needs to be performed manually. We are currently working on an R package that will do both the prediction of expression levels and the association with the phenotype.

All the scripts used to develop the prediction models can be found here:

https://github.com/hwheeler01/PrediXmod





Supported operating systems:
Linux or Mac Os.

The following arguments are allowed, with default values as follows

1. genelist: list of genes. By default it will use all available genes in model database
2. dosages: imputed genotype file path. Default value: 'data/dosages/'
3. dosage_prefix: prefix of dosage file. Default value: 'chr' 
4. weights: full name of database. Default value: 'weight.db'
5. output: output file name 'output'

dosage file format:
- columns are snpid rsid position allele1 allele2 MAF dosage_1 ..... dosage_n 
- dosage for each person refers to the number of alleles for the 2nd allele listed (between 0 and 2)
- it is expected that there will be one file per chromosome

USAGE:
> ./predict_gene_expression.py  --dosages dosagefile_path  --dosages_prefix chr --weights Skin-SunExposed\(Lowerleg\)_0.5.db --output output

** to speed up the process the dosage files can be filtered to SNPs in HapMapSnpsCEU.list.gz.

hapmapSnpsCEU.list.gz: List of SNPs used to develop the models

predict_gene_expression.py: Python script to compute predicted expression
