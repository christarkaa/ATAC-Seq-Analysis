#!/usr/bin/env zsh

# Make a new directory for the project and change to the directory
mkdir ATAC-Seq && cd "$_"

# Download datasets
fastq-dump --split-files SRR24135553 SRR24135554 SRR24135555 SRR24135556

# Make directories for quality control, Mapping and Peak Calling
mkdir QC_Reports Mapping Peaks

# Quality control using fastqc
SRR24135553_1.fastq SRR24135553_2.fastq -o QC_Reports
SRR24135554_1.fastq SRR24135554_2.fastq -o QC_Reports
SRR24135555_1.fastq SRR24135555_2.fastq -o QC_Reports
SRR24135556_1.fastq SRR24135556_2.fastq -o QC_Reports

# Summarizing the QC results
multiqc QC_Reports

