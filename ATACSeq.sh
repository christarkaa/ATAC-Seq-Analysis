#!/usr/bin/env zsh

# Make a new directory for the project and change to the directory
mkdir ATAC-Seq && cd "$_"

# Download datasets
fastq-dump --split-files SRR24135553 SRR24135554 SRR24135555 SRR24135556

# Make directories for quality control, Mapping and Peak Calling
mkdir QC_Reports Mapping Peaks

# Quality control using fastqc
fastqc SRR24135553_1.fastq SRR24135553_2.fastq -o QC_Reports
fastqc SRR24135554_1.fastq SRR24135554_2.fastq -o QC_Reports
fastqc SRR24135555_1.fastq SRR24135555_2.fastq -o QC_Reports
fastqc SRR24135556_1.fastq SRR24135556_2.fastq -o QC_Reports

# Summarizing the QC results
multiqc QC_Reports

# Remove Nextera transposase adapter sequences using cudaput
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
-o trimmed_SRR24135553_1.fastq -p trimmed_SRR24135553_2.fastq \
--report=full SRR24135553_1.fastq SRR24135553_2.fastq > SRR24135553_report.txt

cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
-o trimmed_SRR24135554_1.fastq -p trimmed_SRR24135554_2.fastq \
--report=full SRR24135554_1.fastq SRR24135554_2.fastq > SRR24135554_report.txt

cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
-o trimmed_SRR24135555_1.fastq -p trimmed_SRR24135555_2.fastq \
--report=full SRR24135555_1.fastq SRR24135555_2.fastq > SRR24135555_report.txt

cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
-o trimmed_SRR24135556_1.fastq -p trimmed_SRR24135556_2.fastq \
--report=full SRR24135556_1.fastq SRR24135556_2.fastq > SRR24135556_report.txt

# Trimming with Sickle
sickle pe -f trimmed_SRR24135553_1.fastq -r trimmed_SRR24135553_2.fastq \
-t sanger -o quality_trimmed_SRR24135553_1.fastq -p quality_trimmed_SRR24135553_2.fastq \
-s singletons_SRR24135553.fastq -q 20 -l 90

sickle pe -f trimmed_SRR24135554_1.fastq -r trimmed_SRR24135554_2.fastq \
-t sanger -o quality_trimmed_SRR24135554_1.fastq -p quality_trimmed_SRR24135554_2.fastq \
-s singletons_SRR24135554.fastq -q 20 -l 90

sickle pe -f trimmed_SRR24135555_1.fastq -r trimmed_SRR24135555_2.fastq \
-t sanger -o quality_trimmed_SRR24135555_1.fastq -p quality_trimmed_SRR24135555_2.fastq \
-s singletons_SRR24135555.fastq -q 20 -l 90

sickle pe -f trimmed_SRR24135556_1.fastq -r trimmed_SRR24135556_2.fastq \
-t sanger -o quality_trimmed_SRR24135556_1.fastq -p quality_trimmed_SRR24135556_2.fastq \
-s singletons_SRR24135556.fastq -q 20 -l 90

# Mapping
## Build reference genome index
bowtie2-build GRCh38.p14.genome.fa hg38

## Map reference genome index with the trimmed sequences
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135553_1.fastq -2 quality_trimmed_SRR24135553_2.fastq -S Mapping/SRR24135553.sam
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135554_1.fastq -2 quality_trimmed_SRR24135554_2.fastq -S Mapping/SRR24135554.sam
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135555_1.fastq -2 quality_trimmed_SRR24135555_2.fastq -S Mapping/SRR24135555.sam
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135556_1.fastq -2 quality_trimmed_SRR24135556_2.fastq -S Mapping/SRR24135556.sam


