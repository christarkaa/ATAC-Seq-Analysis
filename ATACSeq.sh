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

# Build reference genome index
bowtie2-build GRCh38.p14.genome.fa hg38

# Map reference genome index with the trimmed sequences
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135553_1.fastq -2 quality_trimmed_SRR24135553_2.fastq -S Mapping/SRR24135553.sam
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135554_1.fastq -2 quality_trimmed_SRR24135554_2.fastq -S Mapping/SRR24135554.sam
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135555_1.fastq -2 quality_trimmed_SRR24135555_2.fastq -S Mapping/SRR24135555.sam
bowtie2 -p 10 -x hg38 -1 quality_trimmed_SRR24135556_1.fastq -2 quality_trimmed_SRR24135556_2.fastq -S Mapping/SRR24135556.sam

# Convert SAM to BAM files using samtools
samtools view -@ 20 -S -b Mapping/SRR24135553.sam > Mapping/SRR24135553.bam
samtools view -@ 20 -S -b Mapping/SRR24135554.sam > Mapping/SRR24135554.bam
samtools view -@ 20 -S -b Mapping/SRR24135555.sam > Mapping/SRR24135555.bam
samtools view -@ 20 -S -b Mapping/SRR24135556.sam > Mapping/SRR24135556.bam

# Elimination of mitochondrial reads using samtools
samtools view -hF 4 Mapping/SRR24135553.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135553_uM.bam
samtools view -hF 4 Mapping/SRR24135554.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135554_uM.bam
samtools view -hF 4 Mapping/SRR24135555.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135555_uM.bam
samtools view -hF 4 Mapping/SRR24135556.bam | grep -vF chrM | samtools view -bS > Mapping/SRR24135556_uM.bam

# Restrict to properly-paired reads only
samtools view -bh -f 3 Mapping/SRR24135553_uM.bam > Mapping/filt_SRR24135553_uM.bam
samtools view -bh -f 3 Mapping/SRR24135554_uM.bam > Mapping/filt_SRR24135554_uM.bam
samtools view -bh -f 3 Mapping/SRR24135555_uM.bam > Mapping/filt_SRR24135555_uM.bam
samtools view -bh -f 3 Mapping/SRR24135556_uM.bam > Mapping/filt_SRR24135556_uM.bam

# Sort the BAM files by coordinate
samtools sort -o Mapping/sorted_SRR24135553_uM.bam Mapping/filt_SRR24135553_uM.bam
samtools sort -o Mapping/sorted_SRR24135554_uM.bam Mapping/filt_SRR24135554_uM.bam
samtools sort -o Mapping/sorted_SRR24135555_uM.bam Mapping/filt_SRR24135555_uM.bam
samtools sort -o Mapping/sorted_SRR24135556_uM.bam Mapping/filt_SRR24135556_uM.bam

# Index the sorted BAM files
samtools index Mapping/sorted_SRR24135553_uM.bam
samtools index Mapping/sorted_SRR24135554_uM.bam
samtools index Mapping/sorted_SRR24135555_uM.bam
samtools index Mapping/sorted_SRR24135556_uM.bam

# Count number of alignments for normalization
count_SRR24135553=$(samtools view -c Mapping/sorted_SRR24135553_uM.bam)
count_SRR24135554=$(samtools view -c Mapping/sorted_SRR24135554_uM.bam)
count_SRR24135555=$(samtools view -c Mapping/sorted_SRR24135555_uM.bam)
count_SRR24135556=$(samtools view -c Mapping/sorted_SRR24135556_uM.bam)

# Print the counts for verification
echo "Alignment counts:"
echo "SRR24135553: $count_SRR24135553"
echo "SRR24135554: $count_SRR24135554"
echo "SRR24135555: $count_SRR24135555"
echo "SRR24135556: $count_SRR24135556"

# Define the lowest read count
lowest_count=5693068

# Calculate subsampling ratios
ratio_SRR24135553=$(awk "BEGIN {print $lowest_count/$count_SRR24135553}")
ratio_SRR24135554=$(awk "BEGIN {print $lowest_count/$count_SRR24135554}")
ratio_SRR24135555=1  # No subsampling needed for SRR24135555
ratio_SRR24135556=$(awk "BEGIN {print $lowest_count/$count_SRR24135556}")

# Print calculated ratios for verification
echo "Subsampling Ratios:"
echo "SRR24135553: $ratio_SRR24135553"
echo "SRR24135554: $ratio_SRR24135554"
echo "SRR24135555: $ratio_SRR24135555"
echo "SRR24135556: $ratio_SRR24135556"

# Subsample using samtools
samtools view -h -b -s $ratio_SRR24135553 Mapping/sorted_SRR24135553_uM.bam > Mapping/subsampled_SRR24135553.bam
samtools view -h -b -s $ratio_SRR24135554 Mapping/sorted_SRR24135554_uM.bam > Mapping/subsampled_SRR24135554.bam
cp Mapping/sorted_SRR24135555_uM.bam Mapping/subsampled_SRR24135555.bam  # No subsampling, just copy
samtools view -h -b -s $ratio_SRR24135556 Mapping/sorted_SRR24135556_uM.bam > Mapping/subsampled_SRR24135556.bam

# Index the subsampled BAM files
samtools index Mapping/subsampled_SRR24135553.bam
samtools index Mapping/subsampled_SRR24135554.bam
samtools index Mapping/subsampled_SRR24135555.bam
samtools index Mapping/subsampled_SRR24135556.bam

# Optional: Count the number of alignments in each subsampled file to verify normalization
echo "Post-subsampling alignment counts:"
samtools view -c Mapping/subsampled_SRR24135553.bam
samtools view -c Mapping/subsampled_SRR24135554.bam
samtools view -c Mapping/subsampled_SRR24135555.bam
samtools view -c Mapping/subsampled_SRR24135556.bam

# Fix mate information
samtools sort -n -o Mapping/namesorted_SRR24135553.bam Mapping/subsampled_SRR24135553.bam
samtools sort -n -o Mapping/namesorted_SRR24135554.bam Mapping/subsampled_SRR24135554.bam
samtools sort -n -o Mapping/namesorted_SRR24135555.bam Mapping/subsampled_SRR24135555.bam
samtools sort -n -o Mapping/namesorted_SRR24135556.bam Mapping/subsampled_SRR24135556.bam

# Fix mate information with samtools fixmate
samtools fixmate -m Mapping/namesorted_SRR24135553.bam Mapping/fixed_SRR24135553.bam
samtools fixmate -m Mapping/namesorted_SRR24135554.bam Mapping/fixed_SRR24135554.bam
samtools fixmate -m Mapping/namesorted_SRR24135555.bam Mapping/fixed_SRR24135555.bam
samtools fixmate -m Mapping/namesorted_SRR24135556.bam Mapping/fixed_SRR24135556.bam

# Sort by coordinates before marking duplicates
samtools sort -o Mapping/sorted_noDup_SRR24135553.bam Mapping/fixed_SRR24135553.bam
samtools sort -o Mapping/sorted_noDup_SRR24135554.bam Mapping/fixed_SRR24135554.bam
samtools sort -o Mapping/sorted_noDup_SRR24135555.bam Mapping/fixed_SRR24135555.bam
samtools sort -o Mapping/sorted_noDup_SRR24135556.bam Mapping/fixed_SRR24135556.bam

# Remove PCR duplicates after sorting by coordinates
samtools markdup -r Mapping/sorted_noDup_SRR24135553.bam Mapping/noDup_SRR24135553.bam
samtools markdup -r Mapping/sorted_noDup_SRR24135554.bam Mapping/noDup_SRR24135554.bam
samtools markdup -r Mapping/sorted_noDup_SRR24135555.bam Mapping/noDup_SRR24135555.bam
samtools markdup -r Mapping/sorted_noDup_SRR24135556.bam Mapping/noDup_SRR24135556.bam

# Index the BAM files after removing duplicates
samtools index Mapping/noDup_SRR24135553.bam
samtools index Mapping/noDup_SRR24135554.bam
samtools index Mapping/noDup_SRR24135555.bam
samtools index Mapping/noDup_SRR24135556.bam

# Convert BAM to BED files using bedtools
bedtools bamtobed -i Mapping/noDup_SRR24135553.bam > Mapping/SRR24135553.bed
bedtools bamtobed -i Mapping/noDup_SRR24135554.bam > Mapping/SRR24135554.bed
bedtools bamtobed -i Mapping/noDup_SRR24135555.bam > Mapping/SRR24135555.bed
bedtools bamtobed -i Mapping/noDup_SRR24135556.bam > Mapping/SRR24135556.bed

# Peak calling macs2
macs2 callpeak -t Mapping/SRR24135553.bed -n SRR24135553 -q 0.05  -f BEDPE -g hs -B --keep-dup all --SPMR --call-summits --outdir Peaks
macs2 callpeak -t Mapping/SRR24135554.bed -n SRR24135554 -q 0.05  -f BEDPE -g hs -B --keep-dup all --SPMR --call-summits --outdir Peaks
macs2 callpeak -t Mapping/SRR24135555.bed -n SRR24135555 -q 0.05  -f BEDPE -g hs -B --keep-dup all --SPMR --call-summits --outdir Peaks
macs2 callpeak -t Mapping/SRR24135556.bed -n SRR24135556 -q 0.05  -f BEDPE -g hs -B --keep-dup all --SPMR --call-summits --outdir Peaks

# Convert narrowPeak files to BED format
for file in Peaks/*.narrowPeak; do
    base=$(basename $file .narrowPeak)
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4}' $file > Peaks/${base}.bed
done

# Format the bed files
awk '$1 !~ /\.[123]$/' Peaks/SRR24135553_peaks.bed > Peaks/SRR24135553_peaks1.bed
awk '$1 !~ /\.[123]$/' Peaks/SRR24135554_peaks.bed > Peaks/SRR24135554_peaks1.bed
awk '$1 !~ /\.[123]$/' Peaks/SRR24135555_peaks.bed > Peaks/SRR24135555_peaks1.bed
awk '$1 !~ /\.[123]$/' Peaks/SRR24135556_peaks.bed > Peaks/SRR24135556_peaks1.bed

# Merge treatment peaks
cat Peaks/SRR24135553_peaks1.bed Peaks/SRR24135554_peaks1.bed | sort -k1,1 -k2,2n > Peaks/treatment_peaks.bed
bedtools merge -i Peaks/treatment_peaks.bed > Peaks/dmPGE2_peaks.bed

# Merge control peaks
cat Peaks/SRR24135555_peaks1.bed Peaks/SRR24135556_peaks1.bed | sort -k1,1 -k2,2n > Peaks/control_peaks.bed
bedtools merge -i Peaks/control_peaks.bed > Peaks/DMSO_peaks.bed

# Peak annotation using homer
annotatePeaks.pl Peaks/dmPGE2_peaks.bed hg38 > Peaks/dmPGE2_peaks.txt
annotatePeaks.pl Peaks/DMSO_peaks.bed hg38 > Peaks/DMSO_peaks.txt

# Convert bedgraph to bigwig files
## Remove the lines in the bedgraph files that are not chromosome names
awk '$1 !~ /\.[123]$/' Peaks/SRR24135553_treat_pileup.bdg > Peak_Calling/filtered.SRR24135553_treat_pileup.bdg
awk '$1 !~ /\.[123]$/' Peaks/SRR24135554_treat_pileup.bdg > Peak_Calling/filtered.SRR24135554_treat_pileup.bdg
awk '$1 !~ /\.[123]$/' Peaks/SRR24135555_treat_pileup.bdg > Peak_Calling/filtered.SRR24135555_treat_pileup.bdg
awk '$1 !~ /\.[123]$/' Peaks/SRR24135556_treat_pileup.bdg > Peak_Calling/filtered.SRR24135556_treat_pileup.bdg

## Get chrom.sizes using samtools
fetchChromSizes hg38 > hg38.chromSizes

# Convert bedGraph to bigWig
bedGraphToBigWig Peaks/filtered.SRR24135553_treat_pileup.bdg hg38.chromSizes Peaks/filtered.SRR24135553.bw
bedGraphToBigWig Peaks/filtered.SRR24135554_treat_pileup.bdg hg38.chromSizes Peaks/filtered.SRR24135554.bw
bedGraphToBigWig Peaks/filtered.SRR24135555_treat_pileup.bdg hg38.chromSizes Peaks/filtered.SRR24135555.bw
bedGraphToBigWig Peaks/filtered.SRR24135556_treat_pileup.bdg hg38.chromSizes Peaks/filtered.SRR24135556.bw
