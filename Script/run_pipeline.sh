#!/usr/bin/env bash
set -euo pipefail

echo "SARS-CoV-2 Viral Genomics Analysis Pipeline"

# 1. Environment Setup
echo "Setting up Bioinformatics Environment..."
sudo apt update && sudo apt upgrade -y
sudo apt install -y fastqc bwa samtools bcftools vcftools
conda install -y -c bioconda freebayes

# 2. Create folders
echo "Creating Folders..."
mkdir -p ../Raw_Data ../Outputs ../Documentation ../Script

# 3. Download data
cd ../Raw_Data
echo "Downloading Raw & Reference Genome Data..."
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O MN908947.fasta
cd ..

# 4. Quality Control
echo "Running FASTQC..."
fastqc ./Raw_Data/ERR5743893_1.fastq.gz ./Raw_Data/ERR5743893_2.fastq.gz --outdir ./Outputs/

# 5. Alignment + BAM processing
echo "Extracting fastq files..."
gunzip -f ./Raw_Data/ERR5743893_1.fastq.gz ./Raw_Data/ERR5743893_2.fastq.gz

cd ./Outputs/
echo "Indexing reference genome with BWA..."
bwa index -p MN908947 ../Raw_Data/MN908947.fasta

echo "Aligning reads with BWA MEM..."
bwa mem MN908947 ../Raw_Data/ERR5743893_1.fastq ../Raw_Data/ERR5743893_2.fastq > ERR5743893.sam

echo "Converting SAM to BAM..."
samtools view -@ 20 -S -b ERR5743893.sam > ERR5743893.bam

echo "Sorting & Indexing BAM file..."
samtools sort -@ 4 -o ERR5743893.sorted.bam ERR5743893.bam
samtools index ERR5743893.sorted.bam

echo "Indexing Reference Genome for variant calling..."
samtools faidx ../Raw_Data/MN908947.fasta

# 6. Variant calling
echo "Calling variants with freebayes..."
freebayes -f ../Raw_Data/MN908947.fasta ERR5743893.sorted.bam > ERR5743893.vcf

echo "SARS-CoV-2 Viral Genomics Analysis Pipeline Completed Successfully!"
