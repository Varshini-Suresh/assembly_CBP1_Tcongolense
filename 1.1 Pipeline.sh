#!/bin/bash

#### SCRIPT TO RUN THE DE NOVO ASSEMBLY PIPELINE OF ONT READS OF T.CONGOLENSE SAMPLES ####

# REQUIRED FILES & DIRECTORIES 
raw_reads="/export/jessie3/2588550s/raw_reads"
ref_genome="/datastore/Genomes/Trypanosoma_congolense/TriTrypDB/release-68/whole_genome/genome.fa"
ref_annotation="/datastore/Genomes/Trypanosoma_congolense/TriTrypDB/release-68/annotation/genes.gff"
Illumina_R1="/export/jessie3/gmh5n/WellcomeCentre/FedericaGiordani/Project3058/T.congolense/trimmed_reads/TcWTcg_S86_L001_R1_001.fastq.gz"
Illumina_R2="/export/jessie3/gmh5n/WellcomeCentre/FedericaGiordani/Project3058/T.congolense/trimmed_reads/TcWTcg_S86_L001_R2_001.fastq.gz"
lineage="euglenozoa_odb10"  
Illumina_reads="/export/jessie3/2588550s/illumina_reads/interleaved_reads.fastq.gz"

# EXTERNAL SCRIPTS
Smartdenovo="/export/jessie3/2588550s/packages/smartdenovo/smartdenovo.pl"
QUAST="/home/2588550S/.local/bin/quast.py"
QC_R="/export/jessie3/2588550s/testing/QC/QC_Plots.R"


# CREATE DIRECTORIES TO STORE THE OUTPUTS
QC_output_raw="/export/jessie3/2588550s/QC/Raw"
processed_reads="/export/jessie3/2588550s/processed_reads"
QC_output_pro="/export/jessie3/2588550s/QC/Processed"
assemblies="/export/jessie3/2588550s/assemblies"
busco_draft_summary="/export/jessie3/2588550s/assemblies/busco_draft_assemblies"
busco_polished_summary="/export/jessie3/2588550s/assemblies/busco_polished_assemblies"
annotation="/export/jessie3/2588550s/annotation"
mkdir -p ${processed_reads}
mkdir -p ${assemblies}
mkdir -p ${busco_draft_summary}
mkdir -p ${busco_polished_summary}
mkdir -p ${annotation}


# Interleave the Illumina read pairs into a single FASTQ file as required by Racon Polishing
seqfu ilv -1 ${Illumina_R1} -2 ${Illumina_R2} > interleave.fastq.gz
tr ' ' ':' < interleave.fastq.gz > ${Illumina_reads} # replace whitespace with ':'
rm interleave.fastq.gz


#### MAIN SCRIPT #####

for sample in  Tco_WTCg_1.7yrs Tco_Drug01_AD6 Tco_Drug03_BF3 Tco_Drug02_BC7 Tco_Drug04_CA9 Tco_Drug01_CE8 

do 

  ### QUALITY CONTROL (QC) OF THE RAW T.CONGOLENSE READS ###
  NanoPlot -t 2 --fastq "${raw_reads}/${sample}.fastq.gz" -o "${QC_output_raw}/${sample}/" --loglength --N50 --raw 
  # Run the QC_PLots R script to recreate the QC plot from NanoPlot's extracted data
  Rscript ${QC_R} --colour "#427d59" --logfile "${QC_output_raw}/${sample}/log.txt" "${QC_output_raw}/${sample}/NanoPlot-data.tsv.gz" "${QC_output_raw}/${sample}/QC_PLOT.png"


  ### READ TRIMMING & FILTERING AND QC ###
  # Creating trimmed reads : Using PoreChop to trim out adapters 
  porechop -i ${raw_reads}/${sample}.fastq.gz --format fastq.gz --no_split -o "${processed_reads}/${sample}_trimmed.fastq.gz"
  # Creating filtered reads : filter out the short reads from the trimmed reads 
  gunzip -c  "${processed_reads}/${sample}_trimmed.fastq.gz" | NanoFilt -l 1000 | gzip > "${processed_reads}/${sample}.fastq.gz"


  ### QUALITY CONTROL (QC) OF THE PROCESSED T.CO READS ###
  NanoPlot -t 2 --fastq "${processed_reads}/${sample}.fastq.gz" -o "${QC_output_pro}/${sample}/" --loglength --N50 --raw 
  # Run the QC_PLots R script to recreate the QC plot from NanoPlot's extracted data
  Rscript ${QC_R}/QC_Plots.R --colour "#6540a1" --logfile "${QC_output_pro}/${sample}/log.txt" "${QC_output_pro}/${sample}/NanoPlot-data.tsv.gz" "${QC_output_pro}/${sample}/QC_PLOT.png"


  ### ASSEMBLE A DRAFT GENOME USING SMARTDENOVO ###
  mkdir -p ${assemblies}/${sample}/
  ${Smartdenovo} -p ${sample} -t 2 -c 1 -J 1000 ${processed_reads}/${sample}.fastq.gz > ${assemblies}/${sample}/${sample}.mak
  make -C ${assemblies}/${sample}/ -f ${assemblies}/${sample}/${sample}.mak 
  SD_assembly="${assemblies}/${sample}/${sample}.dmo.cns.fasta"
  cp ${assemblies}/${sample}/${sample}.dmo.cns ${SD_assembly}
  
  
  ### RUN QC ASSESSMENT ON THE DRAFT ASSEMBLY ###
  # Running the QUAST Assessment 
  python ${QUAST} -t 2 -r ${ref_genome} -g ${ref_annotation} -o "${assemblies}/${sample}/QUAST/" ${SD_assembly}
  # Running BUSCO Assessment 
  busco -m genome -i ${SD_assembly} -o "assemblies/${sample}/${sample}_busco_draft" -l ${lineage} 
  # Copy the BUSCO summaries into a single directory for plotting 
	cp ${assemblies}/${sample}/${sample}_busco_draft/short_summary.*.${lineage}.*.txt ${busco_draft_summary}
  
  
  ### POLISHING THE ASSEMBLY WITH ILLUMINA READS ###
  # Aligning the Smartdenovo draft assembly with Illumina reads (interleaved fastq file)
  minimap2 -a ${SD_assembly} ${Illumina_reads} > "${assemblies}/${sample}/${sample}.sam"
  
  # Run Racon Polishing
  polished_assembly="${assemblies}/${sample}/${sample}_polished_assembly.fasta"  
  racon ${Illumina_reads} ${assemblies}/${sample}/${sample}.sam ${SD_assembly} > ${polished_assembly}
  
  
  ### RUN QC ASSESSMENT ON THE POLISHED ASSEMBLY ###
  # Running the QUAST Assessment 
  python ${QUAST} -t 2 -r ${ref_genome} -g ${ref_annotation} -o "${assemblies}/${sample}/QUAST_polished/" ${polished_assembly}
  # Running BUSCO Assessment 
  busco -m genome -i ${polished_assembly} -o "assemblies/${sample}/${sample}_busco_polished" -l ${lineage} 
  # Copy the BUSCO summaries into a single directory for plotting 
	cp ${assemblies}/${sample}/${sample}_busco_polished/short_summary.*.${lineage}.*.txt ${busco_polished_summary}
  
  
  ## Remove any unwanted files ##
  rm ${processed_reads}/${sample}_trimmed.fastq.gz
  rm busco*log 
  
  # Copy the assembly files to a directory for subsequent annotation 
  cp ${assemblies}/${sample}/${sample}_polished_assembly.fasta ${annotation}/${sample}.fasta
  
done 

# Generate BUSCO summary plots 
generate_plot.py -wd ${busco_polished_summary} 
rm busco*log 
