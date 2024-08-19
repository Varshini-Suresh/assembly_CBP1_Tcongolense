# READ LENGTH & QUALITY PLOTS # 

# Import the required packages 
if (!require("argparse")) install.packages("argparse", repos = "http://cran.us.r-project.org")

library("ggplot2")
library(argparse)

parser <- ArgumentParser(description='Plot reads length vs quality from NanoPlot QC results')
parser$add_argument('tsv_file', help='Input the extracted data (tsv format) generated from NanoPlot when using the --raw command')
parser$add_argument('--colour', type='character', help="Specify a colour to be used for plotting under ('') ")
parser$add_argument('output_dir', help="Specify a directory and filename to store the output (.png)")
parser$add_argument('--logfile', help="Specify a filename to store the read length statistics")

args <- parser$parse_args()

#parser$print_help() # To test the help message 

### Function to plot the graph 
plot_qc <- function(tsv, colour, output, logfile){
  # Read in the TSV file as a table 
  read_data = read.table(tsv, header=TRUE)
  
  # Get the min and max read lengths from the file
  readLengths = paste0("Min length: ", min(read_data$length), "\nMax length: ", max(read_data$length))  
  
  # Plot the data 
  plot = ggplot(data=read_data, aes(x=log10(lengths), y=quals)) + 
    geom_point(colour = colour, alpha = 0.3) + 
    ylim(0, 20) + 
    labs(x = "Read length ", y = "Phred quality score", 
         title = "Read length vs quality") + 
    theme(axis.title = element_text(face = "bold", size = 14), 
          title = element_text(face = "bold", size = 16)) 
  
  # Save the plot to a png in the specified directory 
  # Dimensions are specified in inches 
  ggsave(output, plot, height=5, width=7) 
  
  # Write the read lengths to a log file 
  writeLines(readLengths, logfile)
  
}

plot_qc(args$tsv_file, args$colour, args$output_dir, args$logfile)
