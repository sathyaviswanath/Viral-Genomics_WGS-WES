# **SARS-CoV-2 Viral Genomics Analysis Pipeline**

## 🎯 Pipeline Overview
This pipeline analyzes COVID-19 sequencing data (**ERR5743893**) against the Wuhan-Hu-1 reference genome (**MN908947.3**). It performs quality control, read alignment, BAM processing, and variant calling to identify SNPs/indels.[file:1][file:9]

**Pipeline Workflow Summary**: FASTQ → FastQC → BWA-MEM → Samtools → Freebayes → VCF

## 🚀 Quick Start
## 1. Clone & Setup Environment

git clone <https://github.com/sathyaviswanath/Viral-Genomics_WGS-WES.git>

cd Viral-Genomics_WGS-WES

sudo apt update & sudo apt install fastqc bwa samtools bcftools vcftools

conda install -c bioconda freebayes

## 2. Download Data & Run Pipeline

Create directories and download data, see [Documentation/Pipeline.md](Documentation/Pipeline.md)

## 3. View Results

**Quality reports**

Open Outputs/.html

**Alignment coverage**

Samtools depth Outputs/ERR5743893.sorted.bam > coverage.txt

**Variants**

Less Outputs/ERR5743893.vcf

## 📁 Directory Structure

- Viral-Genomics_WGS-WES
  - `README.md` – Quickstart guide
  - `Documentation/` – Detailed docs  
    - `1_Project_Overview.md`  
    - `2_Pipeline.md`  
  - `Raw_Data/` – FASTQ + reference  
  - `Outputs/` – Results (FastQC, BAM, VCF)  
  - `Script/` – Automation scripts


## 🛠️ Technical Pipeline
| Component   | Tool      | Input                      | Output              |
|-------------|-----------|----------------------------|---------------------|
| **QC**      | FastQC    | FASTQ.gz (paired-end)      | HTML reports  |
| **Alignment** | BWA-MEM | FASTQ + Reference FASTA    | SAM file    |
| **Processing** | Samtools | SAM/BAM                   | Sorted/indexed BAM  |
| **Variants** | Freebayes | BAM + Reference           | VCF 4.2 format |

## 📖 Full Documentation
See [Documentation/Pipeline.md](Documentation/Pipeline.md) for:
- Tool explanations (FastQC, BWA, Samtools, BCFtools, Freebayes)
- Step-by-step commands

## 📚 References
- [ENA: ERR5743893](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/)
- [SARS-CoV-2 Reference](https://www.ebi.ac.uk/ena/browser/view/MN908947.3)


