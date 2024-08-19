## FILTERING THE CBP REGION ANNOTATION FILES ##

# Install the required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rtracklayer")

library(rtracklayer)


# Function to filter a GFF file ####
#' @param filename Annotation of the CBP contig (gff3 format) 
#'
#' @return A GFF file filtered to only include annotations where: 
#' (1) the feature type is "protein match"
#' (2) the signature description is "Serine carboxypeptidase"
#' 
#' @return The number of CBP genes found in the contig
#' 

filter_gff = function (gff_file){
  features = readGFF(gff_file)
  filtered_features = subset(features, features$type == "protein_match" & 
                      features$signature_desc =="Serine carboxypeptidase")
  export(filtered_features, paste0("filtered_", gff_file), format = "GFF3")
  return (paste0('Number of CBP genes in ', gff_file, ' are ',nrow(filtered_features)))
}


# Run the function on all the annotation files from each sample 
filter_gff('WT_utg64.gff3') # 13
filter_gff('AD6_utg299.gff3') # 15
filter_gff('BF3_utg2430.gff3') # 2 
filter_gff('CE8_utg235.gff3') # 19 
filter_gff('BC7_utg639.gff3') # 15
filter_gff('CA9_utg763.gff3') # 3 
