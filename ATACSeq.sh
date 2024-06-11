#!/usr/bin/env zsh

# Make a new directory for the project and change to the directory
mkdir ATAC-Seq && cd "$_"

# Download datasets
fastq-dump --split-files SRR24135553 SRR24135555

# Make directories for quality control, Mapping and Peak Calling
mkdir QC_Reports Mapping Peaks

# Quality control using fastqc

