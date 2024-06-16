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

# Convert sam to bam files using samtools
samtools view -@ 20 -S -b Mapping/SRR24135553.sam > Mapping/SRR24135553.bam
samtools view -@ 20 -S -b Mapping/SRR24135554.sam > Mapping/SRR24135554.bam
samtools view -@ 20 -S -b Mapping/SRR24135555.sam > Mapping/SRR24135555.bam
samtools view -@ 20 -S -b Mapping/SRR24135556.sam > Mapping/SRR24135556.bam

# Elimination of mitochondria reads using samtools
samtools view -hF 4 Mapping/SRR24135553.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135553_uM.bam
samtools view -hF 4 Mapping/SRR24135554.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135554_uM.bam
samtools view -hF 4 Mapping/SRR24135555.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135555_uM.bam
samtools view -hF 4 Mapping/SRR24135556.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135556_uM.bam

# Sort the bam files using samtools
samtools sort Mapping/SRR24135553_uM.bam > Mapping/sorted_SRR24135553_uM.bam
samtools sort Mapping/SRR24135554_uM.bam > Mapping/sorted_SRR24135554_uM.bam
samtools sort Mapping/SRR24135555_uM.bam > Mapping/sorted_SRR24135555_uM.bam
samtools sort Mapping/SRR24135556_uM.bam > Mapping/sorted_SRR24135556_uM.bam

# Index the sorted bam files using samtools
samtools index Mapping/sorted_SRR24135553_uM.bam
samtools index Mapping/sorted_SRR24135554_uM.bam
samtools index Mapping/sorted_SRR24135555_uM.bam
samtools index Mapping/sorted_SRR24135556_uM.bam

# Count number of alignments
samtools view -c Mapping/sorted_SRR24135553_uM.bam
samtools view -c Mapping/sorted_SRR24135554_uM.bam
samtools view -c Mapping/sorted_SRR24135555_uM.bam
samtools view -c Mapping/sorted_SRR24135556_uM.bam

# Normalisation
## The number of alignment count for each sample
# 9762071
# 9923188
# 6770216
# 6886528

# Define the lowest read count
lowest_count=6770216

# Calculate subsampling ratios for each sample
ratio_SRR24135553=$(awk "BEGIN {print $lowest_count/9762071}")
ratio_SRR24135554=$(awk "BEGIN {print $lowest_count/9923188}")
ratio_SRR24135555=1  # No subsampling needed for SRR24135555, already at target size
ratio_SRR24135556=$(awk "BEGIN {print $lowest_count/6886528}")

# Print calculated ratios for verification
echo "Subsampling Ratios:"
echo "SRR24135553: $ratio_SRR24135553"
echo "SRR24135554: $ratio_SRR24135554"
echo "SRR24135555: $ratio_SRR24135555"
echo "SRR24135556: $ratio_SRR24135556"

# Subsample using samtools
samtools view -h -b -s 660$ratio_SRR24135553 Mapping/sorted_SRR24135553_uM.bam > Mapping/subsampled_SRR24135553.bam
samtools view -h -b -s 660$ratio_SRR24135554 Mapping/sorted_SRR24135554_uM.bam > Mapping/subsampled_SRR24135554.bam
samtools view -h -b Mapping/sorted_SRR24135555_uM.bam > Mapping/subsampled_SRR24135555.bam  # No subsampling, just copy
samtools view -h -b -s 660$ratio_SRR24135556 Mapping/sorted_SRR24135556_uM.bam > Mapping/subsampled_SRR24135556.bam

# Index the subsampled BAM files
samtools index Mapping/subsampled_SRR24135553.bam
samtools index Mapping/subsampled_SRR24135554.bam
samtools index Mapping/subsampled_SRR24135555.bam  # Index the copied BAM file
samtools index Mapping/subsampled_SRR24135556.bam

# Optional: Count the number of alignments in each subsampled file to verify normalization
echo "Post-subsampling alignment counts:"
samtools view -c Mapping/subsampled_SRR24135553.bam
samtools view -c Mapping/subsampled_SRR24135554.bam
samtools view -c Mapping/subsampled_SRR24135555.bam
samtools view -c Mapping/subsampled_SRR24135556.bam

# Convert bam to bed files
bedtools bamtobed -i Mapping/subsampled_SRR24135553.bam > Mapping/SRR24135553.bed
bedtools bamtobed -i Mapping/subsampled_SRR24135554.bam > Mapping/SRR24135554.bed
bedtools bamtobed -i Mapping/subsampled_SRR24135555.bam > Mapping/SRR24135555.bed
bedtools bamtobed -i Mapping/subsampled_SRR24135556.bam > Mapping/SRR24135556.bed


# Peak calling macs2
macs2 callpeak -t Mapping/SRR24135553.bed -n SRR24135553 --qvalue 0.001 -f BED -g hs --keep-dup all --SPMR --call-summits --outdir Peaks/SRR24135553
macs2 callpeak -t Mapping/SRR24135554.bed -n SRR24135554 --qvalue 0.001 -f BED -g hs --keep-dup all --SPMR --call-summits --outdir Peaks/SRR24135554
macs2 callpeak -t Mapping/SRR24135555.bed -n SRR24135555 --qvalue 0.001 -f BED -g hs --keep-dup all --SPMR --call-summits --outdir Peaks/SRR24135555
macs2 callpeak -t Mapping/SRR24135556.bed -n SRR24135556 --qvalue 0.001 -f BED -g hs --keep-dup all --SPMR --call-summits --outdir Peaks/SRR24135556
