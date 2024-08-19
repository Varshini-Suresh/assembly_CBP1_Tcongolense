#### SCRIPT TO EXTRACT CONTIG CONTAINING THE CBP REGION FROM THE GENOME ASSEMBLY ####

## Install the required packages ####
if (!require('devtools')) install.packages('devtools'); library('devtools')
install_github("GrahamHamilton/pipelineTools")
library("pipelineTools", "parallel")

# Run the run_blast function from pipelineTools, located in another directory   
source("/export/jessie3/2588550s/testing/ANNOTATION/run_blast.R")


# Links for BLAST path  ####
path_db <- "/software/blast-v2.13.0/bin/makeblastdb"
blast.path <- "/software/blast-v2.13.0/bin/blastn"

# Link to query CBP1 sequence from the reference T.congolense genome ####
query <- "/export/jessie3/gmh5n/WellcomeCentre/FedericaGiordani/NanoporeData/cbp.fa"

# Set the sample name (cell line) here, eg. 'AD6' ####
sample <- "WT"
  

# Make a custom BLAST database ####
blast_db_cmds <- run_blast(input = paste0("/export/jessie3/2588550s/annotation/",sample,".fasta"),
                           title = sample,
                           database = "nucl",
                           parallel = TRUE,
                           cores = 6,
                           execute = TRUE,
                           blast = path_db)

# BLAST search the assembly with the query CBP1 sequence ####
blast_cmds <- run_blast(input = query,
                        output = paste0(sample,".blast"),
                        database = paste0(sample,".fasta"),
                        format = 6,
                        parallel = TRUE,
                        cores = 6,
                        execute = TRUE,
                        blast = blast.path)

# Remove any intermediate files 
intermediate <- list.files(pattern = paste0(sample,".fasta.n"))
file.remove(intermediate)

# This will generate a DB containing the contig that matched the query CBP1 sequence

# Next is to run the blast parser script to generate the contig sequence 
# python BLAST_parser.py -b ${sample}.blast -a {sample}.fasta


