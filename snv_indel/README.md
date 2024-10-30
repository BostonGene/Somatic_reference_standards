## Description of content

## ID_snv_indel.tsv
Reference set of SNV and INDEL in .TSV format. The dataset contains the following columns: 
*Hugo_Symbol* - the gene symbol acciording to Gencode43
*Chromosome, Start_Position, End_Position* - genome coordinates according to hg38
*Reference_Allele, Tumor_Seq_Allele2* - reference and alternative allele in standard .MAF format
*Variant_Classification* - annotation according to incluence on a protein level
*Variant_Type* - classification into SNP/DEL/INS
*Tumor_VAF_median* - variant allele frequency expected to be at 100% tumor cell line
*VAF_STDEV* - standard deviation of variant callele frequency, can be used to additionally verify the reference set

## ID__non_mutated_genes.tsv
Regions of the genes without SNV/indel mutations. Can be used as a set for precision assessment.The dataset contains the following columns: 
*Chromosome, Start_Position, End_Position, Region* - genome coordinates according to hg38
*Strand* - forward or reverse
*Hugo_Symbol* - the gene symbol acciording to Gencode43
*Type* - type of a region - exon/intron/CDS
*Number* - number of exon
*Ensemble, RefSeq* - ID of a transcript
*QC content, Mapping Quallity* - features of a sequence, can be used to narrow the tested regions 
