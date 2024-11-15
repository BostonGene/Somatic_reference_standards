## For SNV/Indel
Use the script snv_indel_metrics.py to compare mutation data set (SNV + Indel) to reference standard. The output is a .txt file containing Sensitivity and Precision assessments. 
The input file should be in a standard .maf format (columns Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2 are required). 

The example of launch:
` python snv_indel_metrics.py --input test_sample.maf --target_bed test_region.bed --reference_sensitivity test_reference_sensitivity.tsv --reference_precision test_reference_precision.tsv --purity 100 --type SNP --output test_output.txt `

Agruments:

` --input, .maf file with annotated variants `
` --target_bed, Input target.bed file `
` --reference_sensitivity, Reference file for sensitivity `
` --reference_precision, Reference file for precision `
` --purity, Purity of the sample for analysis, must be eiter 10, 20, 30, 50, 75, 100 % `
` --type, Type of the events for analysis `
` --output, Output .txt file with performance characteristics `

## For CNV
Use the script cnv_metrics.py to compare normalized reference statuses of genes for either of  5 model cell lines to the output of CNV calling pipeline. The output is a .txt file containing Sensitivity and Specificity assassments.

The example of launch:
` python cnv_metrics.py --input test-cna-normalized.txt --reference test_reference.tsv --genes test_genes.txt --cell_line NCI-H1770 --output test_output.txt `

Arguments:
` --input, Input .tsv file with normalized variants `
` --reference, Reference file for cnv and neutral genes `
` --genes, Target genes for testing `
` --cell_line, Cell line for testing [COLO829, HCC1143, HCC1395, HCC1937, NCI-H1770] `
` --output, Output .txt file with performance characteristics `
