#!/bin/bash

# Script to extract the CBP1 gene sequence from the GFF coordinates 
# Outputs: 
# (1) FASTA file containing only the CBP1 sequences and 
# (2) A log file showing the lengths of each CBP1 copy present

# KEY
# List of samples and their respective contigs containing the CBP region
# Tco_NoDrug_WT- utg64
# Tco_Drug01_AD6 - utg299
# Tco_Drug02_BF3 - utg2430
# Tco_Drug03_BC7 - utg639
# Tco_Drug04_CA9 - utg763
# Tco_Drug01_CE8 - utg235

# RELEVANT DIRECTORIES
annotation="/export/jessie3/2588550s/annotation"
GFF="/export/jessie3/2588550s/annotation/GFF"

# CREATE OUTPUT DIRECTORY
CBP_out_dir="/export/jessie3/2588550s/annotation/CBP_sequences"
mkdir -p ${CBP_out_dir}

# SET THE VARIABLES 
# Change the sample and contig as required. 
sample="BC7"
contig="utg2430"


# Set the output filenames 
CBP_fasta="${CBP_out_dir}/${sample}_CBP.fasta"
logfile="${CBP_out_dir}/${sample}_${contig}_log.txt"
filt_CBP_fasta="${CBP_out_dir}/filtered_${sample}_CBP.fasta"


## MAIN SCRIPT ##
# Extract sequence of each CBP gene copies from the genomic coordinates (GFF file)
bedtools getfasta -fi ${annotation}/${contig}.fasta -fo ${CBP_fasta} -bed ${GFF}/filtered_${sample}_${contig}.gff3

# Print the number of CBP gene copies in the sample 
copies=`grep ">" ${CBP_fasta} | wc -l` 

# Print message on std on running the script
echo "The ${sample} cell line has ${copies} copies of CBP1 genes in the ${contig} contig." 


# FILTER OUT SHORTER FRAGMENTS # 
# Remove short sequences (< 1000) as they are only fragments and not full CBP1 genes 
bioawk -c fastx '{ if(length($seq) > 1000) { print ">"$name; print $seq }}'> ${filt_CBP_fasta} ${CBP_fasta}

# Print number of genes present after filtering out fragments
copies=`grep ">" ${filt_CBP_fasta} | wc -l`

# Use bioawk package to count the sequence length and write it to a logfile 
bioawk -c fastx '{print $name, length($seq)}' ${filt_CBP_fasta} > ${logfile}

# Print message on std on running the script
echo "The ${sample} cell line now has ${copies} copies of CBP1 genes in the ${contig} contig after filtering out shorter fragments (<1000 bases)." 

####
