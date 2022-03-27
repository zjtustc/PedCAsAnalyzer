# PedCAsAnalyzer
PedCAsAnalyzer is a helpful method for CA(chromosomal aneuploidy) diagnosis and CA type analysis by analyzing mutations from WES data. It fills the gap between the genetic diagnosis of CA and monogenic variants, which improves the disease diagnostic yield and efficiency

## input
**1** VCF files
* mode 1: If only proband VCF file provided, CA diagnosis function is avaiable
* mode 2: If trios VCF file provided, CA diagnosis and CA type analysis function is avaiable

**2** PED file
* PED file is used to describe the sample relationship, details in [PED](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).   
**Notice:** sample ID in PED file should be consistent with VCF name

## Command line usage
* -task task id
* -ped PED file
* -chr chromosomes analysis, for example: chr18 chr21 default:chrX
* -plot Y----plot, N----don't plot
* -MAF parent sepcific marker MAF threshold, default='0.02'

* example. To analysis chr21 and chrX aneuploidy for Pedigree named PK662
  python PedCAsAnalyzer.py -task PK662 -ped PK662.ped -chr chr21
  
