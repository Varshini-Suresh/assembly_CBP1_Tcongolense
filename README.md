# Assembly and Investigation of the Serine Carboxypeptidase (_CBP1_) gene cluster in _Trypanosoma congolense_
## Project Background & Aims
This study aimed to use Oxford Nanopore Technology (ONT) to construct a _de novo_ genome assembly of the wild-type (WT) and five resistant T. congolense cell lines and compare the _CBP1_ genes between the cell lines to gather a better understanding of _CBP1_ cluster region deletions and resistance mechanism of _T. congolense_ against the antitrypanosomal prodrugs.


## A set of scripts written for de novo assembly and investigation of _CBP1_ gene clusters in _T. congolense_ 
**1.1	Pipeline.sh** <br/>
The Bash Script of the entire de novo assembly pipeline done on a for loop on each sample in the study. For the time being, the sample read files (FASTQ) were hard-coded, but the script can be changed to take in the sample ONT reads (FASTQ file format) as command-line options from the user. In this case, the sample read filenames were changed to remove identity of the drug compounds used for each resistant _T. congolense_ cell line. <br/>

**1.2	QC_Plots.R**<br/>
The custom R script was written to use the extracted data from NanoPlot results to make plots without subsetting the data. The script also uses the extracted data to determine the lengths of the shortest and longest reads from the given FASTQ samples. The command-line usage is:<br/>
```Rscript QC_Plots.R [-h] [--colour COLOUR] [--logfile LOGFILE] tsv_file output_dir```<br/>

**1.3	CBP_contig_extraction.R**<br/>
This R Script was written in conjunction to the ‘run_blast’ function from the ‘run_blast.R’ script (https://github.com/GrahamHamilton/pipelineTools/blob/master/R/run_blast.R) available on the pipelineTools package GitHub page.  The script creates a BLAST database using the whole assembly and then reruns the same ‘run_blast’ function with the query CBP1 sequence to extract the contig from the assembly that matches the query sequence.<br/>

**1.4	GFF_filter.R**<br/>
A custom function called “filter_gff” was written on R using the ‘rtracklayer’ package to convert the annotation file (GFF) into a dataframe, which was filtered to only include the _CBP1_ gene annotations in the GFF file. The function also calculates the number of sequences in the contig that were annotated as CBP1 genes.<br/>

**1.5	Extract_CBP1.sh**<br/>
This Bash Script was written to extract the _CBP1_ gene sequences from the entire contig sequence files and save it as a FASTA file. The script also filters out any shorter fragments (<1000 bases) in the CBP1 cluster that also happens to be annotated as a _CBP1_ gene. The lengths of each _CBP1_ gene from the FASTA file is printed out as a separate log file.<br/>

**1.6	MSA_phylo_tree.R**<br/>
This R script was written using several packages to conduct the multiple sequence alignment between the different _CBP1_ sequences. The script starts by merging the input FASTA files from multiple samples (cell lines) and then conducts an MSA. From the pairwise matrices, the phylogenetic trees and plotted and saved in the working directory.<br/>

