#### MULTIPLE SEQUENCE ALIGNMENT AND PHYLOGENETIC TREE - CBP1 REGION ####

## Install the required packages ####
# BiocManager packages
biopackages <- c("msa", "Biostrings", "pwalign", "seqinr", "treeio")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(biopackages)
lapply(biopackages, library, character.only = TRUE)

# Github packages 
devtools::install_github("YuLab-SMU/ggmsa")
library("ggmsa")

# Other packages
pkgs <- c("dplyr", "tidyr", "xfun", "ape", "ggtree", "tinytex", "stats")
if (!requireNamespace(pkgs, quietly = TRUE))
  install.packages(pkgs)
lapply(pkgs, library, character.only = TRUE)


# Set the sample names, files & directories ####
# Set the working directory
setwd("/export/jessie3/2588550s/annotation/CBP_sequences/MSA")

# Enter the samples to be compared in sample1_sample2 format: 
# Example "WT_vivax"
comparison="pubWT_BC7"

# Set the sample fasta file directories" 
sample_dir=(paste0("./",comparison))
merged_fasta = paste0(sample_dir, "/merged_", comparison, ".fasta")
msa_file=paste0(sample_dir,"/",comparison,".pdf")


# Make a sample directory on command line and copy the respective FASTA files onto the directory
# Example: 
# mkdir ${comparison}
# cp ${sample_1}.fasta ${sample_2}.fasta ${comparison}

# Merge the fasta files and write out the file 
files <- list.files(sample_dir, pattern = "_CBP.fasta", recursive = TRUE, full.names = TRUE)
merged_file <- read_all(files)
write.table(merged_file, file=merged_fasta, col.names = FALSE, row.names = FALSE, quote = FALSE)


### Run the Multiple Sequence Analysis ####
# Convert into a XStringSet object for MSA 
fasta <- readDNAStringSet(merged_fasta) 

# Run the MSA in default settings (ClustalW)
multiple_alignment <- msa(fasta)

# Convert MSA to be used by other packages 
multiple_alignment_converted <- msaConvert(multiple_alignment, type="seqinr::alignment")


### Build the Phylogenetic Tree ####
# Get distances
d <- dist.alignment(multiple_alignment_converted, "identity")
tree <- njs(d)

# Set the "group" label, by splitting the FASTA header by "_" to get sample group  
group <- as.treedata(tree) %>%
  as_tibble() %>%
  dplyr::select(label) %>%
  separate_wider_delim(label, names = c("Group", NA), delim = "_", cols_remove = FALSE)

# Create the treeio data 
treeio_data <- as.treedata(tree) %>%
  as_data_frame() %>%
  as.treedata() %>%
  full_join(group, by = 'label')
options(ignore.negative.edge=TRUE)

# Plot the treeio data using ggtree 
ggplot(treeio_data, aes(x, y)) + 
  geom_tree() + 
  theme_tree() + 
  geom_tippoint(aes(color=group$Group), size = 3) +
  geom_tiplab(align=F, hjust = -0.1) +
  scale_x_continuous(expand = expansion(mult = 0.5))

# Save the plot 
ggsave(paste0(sample_dir,"/","filtered_",comparison,"_tree.png"), height = 10, width = 8)


### Plot the MSA for the sequences ####

# Export the MSA into a PDF format 
msaPrettyPrint(multiple_alignment,
               file = msa_file,
               output = "pdf",
               showNames = "left",
               showNumbering = "right",
               showLogo = "top",
               showConsensus = "bottom",
               logoColors = "rasmol",
               verbose = FALSE,
               askForOverwrite = FALSE)


