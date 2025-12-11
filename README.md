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

Create directories and download data, see [Documentation/Pipeline.md](Documentation/Pipeline.md) & [Script](Script/run_pipeline.sh)

## 3. View Results

**Quality reports**

Open Outputs/.html

**Alignment coverage**

Samtools depth Outputs/ERR5743893.sorted.bam > coverage.txt

**Variants**

Less Outputs/ERR5743893.vcf

## 🛠️ Technical Pipeline and Tools Explanation
| Component   | Tool      | Input                      | Output              |
|-------------|-----------|----------------------------|---------------------|
| **QC**      | FastQC    | FASTQ.gz (paired-end)      | HTML reports  |
| **Alignment** | BWA-MEM | FASTQ + Reference FASTA    | SAM file    |
| **Processing** | Samtools | SAM/BAM                   | Sorted/indexed BAM  |
| **Variants** | Freebayes | BAM + Reference           | VCF 4.2 format |

## Tool Explanations
### 1. FastQC

**Purpose:** Generates quality control reports for raw sequencing reads, assessing per-base quality, GC content, sequence duplication, and adapter contamination.

- Essential first step to identify data issues before downstream analysis

- Produces HTML reports with summary statistics and graphs​

### 2. BWA (Burrows-Wheeler Aligner)

**Purpose:** Aligns short sequencing reads to a reference genome using the BWA-MEM algorithm optimized for Illumina paired-end reads.

- Handles mismatches, gaps, and complex mapping scenarios efficiently

- Outputs SAM format with alignment coordinates and mapping quality scores​

### 3. Samtools(Software Asset Management Tools)

**Purpose:** Suite of utilities for manipulating SAM/BAM files and reference genomes.

- Converts formats (SAM↔BAM), sorts alignments, indexes files for fast access

- faidx creates random access indexes for FASTA files used in variant calling​

### 4. Freebayes

**Purpose:** Haplotype-based variant caller for discovering SNPs and indels from BAM alignments.

- Population genetics model considers allele frequencies and mapping quality

- Installed via Bioconda for compatibility with bioinformatics environments.

### 5. BCFtools

**Purpose:** High-performance toolkit for manipulation and analysis of VCF/BCF variant files.

- Filters variants by quality, depth, allele frequency (bcftools view -i 'QUAL>20')

- Merges multiple VCFs, generates summaries, and performs statistical analysis

- Essential for post-variant calling processing and quality control of Freebayes output

## 📖 Documentation
See [Documentation/Pipeline.md](Documentation/Pipeline.md) for:
- Step-by-step commands

See [Documentation/Troubleshooting.md](Documentation/3.Troubleshooting.md) for:
- Common issues

## 📚 References
- [ENA: ERR5743893](https://www.ebi.ac.uk/ena/browser/view/ERR5743893)
- [SARS-CoV-2 Reference](https://www.ebi.ac.uk/ena/browser/view/MN908947.3)


