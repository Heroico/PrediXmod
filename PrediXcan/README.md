PrediXcan
=========

PrediXcan is a gene-based association test that prioritizes genes that are likely to be causal for the phenotype. 

Reference
---------

Under revision

Instructions
------------

Input: 
- genotype (user provided)
- phenotype (user provided)
- transcriptome prediction model (sqlite db)
 	tissue: Whole Blood (default)
	source: DGN (default)
	model: Elastic Net 0.50 (default)

Below is a link to the files needed to run PrediXcan. The folder contains the prediction model for sun exposed skin expression levels (SQLite database), a python script to compute predicted expression levels, and the cross validated performance measures for each gene. 

https://app.box.com/s/uzmgllg2pba90ame3kiwd6nrdglgnsp8

For now the script only computes the predicted expression levels. The association with the expression levels must be done manually. We are developing an R package that we will make available when ready.

All the scripts used to develop the prediction models can be found here:

https://github.com/hwheeler01/PrediXmod

Please let me know if there is anything else you need.

