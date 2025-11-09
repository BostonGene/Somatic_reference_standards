## Description of content

## cell_lines_reference_CNV.whole_exome.tsv
Reference set of Copy Number Variantions (amplification and deletions) in .TSV format. The dataset contains the following columns:
*Gene* - the gene symbol acciording to Gencode43
*HCC1143, COLO829, HCC1937, HCC1395, NCI-H1770* - ID of cell lines in dataset

Normalized copynumbers are provided for each gene for each cell line. The dataset contains:
*0* - genes with copynumber in tumor sample equal to ploidy of the tumor genome,
*-1/+1* - genes with sligntly altered copynumer state compared to ploidy,
*+2* - high-level amplifications, generally, with copynumbers greater then two times ploidy of a tumor sample,
*-2* - gene loss in a tumor sample
