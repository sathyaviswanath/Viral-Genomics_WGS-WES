#!/usr/bin/env bash
set -e

# 1. Environment Setup
sudo apt update & sudo apt upgrade
sudo apt install fastqc bwa samtools bcftools vcftools
conda install -c bioconda freebayes

# 2. Create folders
mkdir -p ../Raw_Data ../Outputs ../Documentation ../Script

# 3. Download data
cd ../Raw_Data
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O MN908947.fasta
cd ..

# 4. QC
fastqc ./Raw_Data/ERR5743893_1.fastq.gz ./Raw_Data/ERR5743893_2.fastq.gz --outdir ./Outputs/

# 5. Alignment + BAM processing
gunzip ./Raw_Data/ERR5743893_1.fastq.gz ./Raw_Data/ERR5743893_2.fastq.gz
cd Outputs/
bwa index -p MN908947 ../Raw_Data/MN908947.fasta
bwa mem MN908947 ../Raw_Data/ERR5743893_1.fastq ../Raw_Data/ERR5743893_2.fastq > ERR5743893.sam
samtools view -@ 20 -S -b ERR5743893.sam > ERR5743893.bam
samtools sort -@ 4 -o ERR5743893.sorted.bam ERR5743893.bam
samtools index ERR5743893.sorted.bam
samtools faidx ../Raw_Data/MN908947.fasta

# 5. Variant calling
freebayes -f ../Raw_Data/MN908947.fasta ERR5743893.sorted.bam > ERR5743893.vcf
