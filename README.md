# GangSTR Pipeline-wdl
Workflows for detecting short tandem repeat expansions using GangSTR (https://github.com/gymreklab/GangSTR) on Terra Platform.These workflows are designed for Whole Genome Sequencing starting from mapped bam/cram files. 

The pipeline includes following tasks:

1.Genotype tandem repeats provided in the catalog using GangSTR
2.Compress the GangSTR output VCF file and ganerate vcf index
3.Filter GangSTR output using dumpSTR (https://trtools.readthedocs.io/en/latest/source/dumpSTR.html)
4.Convert filtered vcf output into compressed tsv file 


