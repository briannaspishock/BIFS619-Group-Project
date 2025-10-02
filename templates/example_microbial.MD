# Group Project: RNA-Seq Workflow (Microbial)

## Preparation
- **GitHub Setup**
  - Team members: Brianna Spishock, Youngju Jeon, Brendan O'Brien 

## Data Retrieval
- **NCBI SRA accessions used:**
  - SRR11998473 -starvation
  - SRR11998467 -oxidative
  - SRR11998457 -acidic


- **Commands used to download data:**
  ```bash

# 1. Download SRA files
prefetch SRR11998457
prefetch SRR11998467
prefetch SRR11998473

# 2. Extract paired-end FASTQ files into condition folders

# Acidic stress
mkdir -p acidic
fasterq-dump SRR11998457 --split-files --threads 8 --outdir acidic

# Oxidative stress
mkdir -p oxidative
fasterq-dump SRR11998467 --split-files --threads 8 --outdir oxidative

# Starvation stress
mkdir -p starvation
fasterq-dump SRR11998473 --split-files --threads 8 --outdir starvation


## Quality Control (QC)
- **Tools Used:** 
  - (Example) FastQC
  - (Example) MultiQC
- **Tasks Performed:**  
  - (Example) Generated per-sample QC reports  
  - (Example) Summarized read quality, GC content, adapter contamination, duplication  
- **Deliverables:**  
  - (Example) Insert QC plots (e.g., per base quality, adapter content)  
  - (Example) Insert Table: Raw read counts, duplication rates  
- **Interpretation:**  
  - Brief narrative of overall quality  
  - Any concerns identified in the QC step  

## Read Cleaning
- **Tools Used:** 
- **Tasks Performed:**  
- **Deliverables:**  
- **Interpretation:**  

