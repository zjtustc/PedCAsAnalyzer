# PedCAsAnalyzer
PedCAsAnalyzer is a helpful method for CA(chromosomal aneuploidy) diagnosis and CA type analysis by analyzing mutations from WES data. It fills the gap between the genetic diagnosis of CA and monogenic variants, which improves the disease diagnostic yield and efficiency

## input
**1** VCF files
* mode 1: If only proband VCF file provided, CA diagnosis function is avaiable
* mode 2: If trios VCF file provided, CA diagnosis and CA type analysis function is avaiable

**2** Ped file
* Ped file is used to describe the sample relationship, details in [Ped]:<https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format>
