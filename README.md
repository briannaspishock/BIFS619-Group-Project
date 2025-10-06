# Group Project: RNA-Seq Workflow (Microbial)

## Preparation
- **GitHub Setup**
  - Team members: Brianna Spishock, Youngju Jeon, Brendan O'Brien 

## Data Retrieval
- **NCBI SRA accessions used:**
  - SRR11998473 -starvation
  - SRR11998467 -oxidative
  - SRR11998457 -acidic

Bioproject PRJNA638918 "RNA Atlas of Bacterial Human Pathogens Uncovers Stress Dynamics"

## Download Salmonella RNA-seq Data

```bash
# Navigate to analysis directory
cd ~/groupproject

# Acidic stress
mkdir -p acidic
prefetch SRR11998457
fasterq-dump SRR11998457 --split-files --threads 8 --outdir acidic

# Oxidative stress
mkdir -p oxidative
prefetch SRR11998467
fasterq-dump SRR11998467 --split-files --threads 8 --outdir oxidative

# Starvation stress
mkdir -p starvation
prefetch SRR11998473
fasterq-dump SRR11998473 --split-files --threads 8 --outdir starvation
```


## Quality Control (QC)
- **Tools Used:** 
  - FastQC
  - MultiQC
 
```bash
# Create an output folder for FastQC reports
mkdir -p qc_reports

# Run FastQC on all FASTQ files inside each condition folder
fastqc acidic/*.fastq -o qc_reports
fastqc oxidative/*.fastq -o qc_reports
fastqc starvation/*.fastq -o qc_reports
```

```bash
# Put all reports together using MultiQC
cd qc_reports
multiqc .
# Generate QC plots 
multiqc . --export
```

- **Tasks Performed:**  
  - Ran FastQC on all paired-end FASTQ files across three conditions (acidic, oxidative, starvation)
  - Generated per-sample QC reports (`*_fastqc.html`) with metrics for sequence quality, GC content, duplication, and adapter contamination  
  - Ran MultiQC to put all FastQC results into a single summary report (`multiqc_report.html`)  

- **Deliverables:**
### 1-2 QC plots raw data
  - [MultiQC Report]()
  - [Adapter Content Plot]()
  - [Sequence Duplication Levels]()
  - [Per Sequence GC Content]()

### Raw Read Counts and Duplication Rates (from multiqc_general_stats.txt)

| Sample        | Condition   | Raw Reads | Duplication % | GC % | Avg Length (bp) |
|---------------|-------------|-----------|--------------|------|------------------|
| SRR11998457_1 | Acidic      | 22.7 M    | 89.3%         | 50%  | 143             |
| SRR11998457_2 | Acidic      | 22.7 M    | 90.7%         | 50%  | 151             |
| SRR11998467_1 | Oxidative   | 10.2 M    | 87.5%         | 51%  | 143             |
| SRR11998467_2 | Oxidative   | 10.2 M    | 88.4%         | 51%  | 151             | 
| SRR11998473_1 | Starvation  | 5.4 M     | 88.1%         | 51%  | 143             | 
| SRR11998473_2 | Starvation  | 5.4 M     | 88.4%         | 50%  | 151             | 

- **Interpretation:**  
  - Read counts: Range from 22.7M in acidic to 5.4M in starvation. This reflects sequencing depth across different conditions.
  - Duplication: 87.5-90.7%
  - GC content: Stable at 50-51%
  - Read length: forward reads are 143bp and reverse reads are 151bp, consistent with expected paired-end design.

## Read Cleaning
- Tools Used:
  - Fastp
  - MutliQC

 ```bash
#Acidic
fastp \
  -i acidic/SRR11998457_1.fastq \
  -I acidic/SRR11998457_2.fastq \
  -o trimmed/SRR11998457_1.clean.fastq.gz \
  -O trimmed/SRR11998457_2.clean.fastq.gz \
  --detect_adapter_for_pe \
  --disable_quality_filtering \
  --length_required 30 \
  --thread 8 \
  --html qc_reports_trimmed/SRR11998457_fastp.html \
  --json qc_reports_trimmed/SRR11998457_fastp.json 

#Oxidative
fastp \
  -i oxidative/SRR11998467_1.fastq \
  -I oxidative/SRR11998467_2.fastq \
  -o trimmed/SRR11998467_1.clean.fastq.gz \
  -O trimmed/SRR11998467_2.clean.fastq.gz \
  --detect_adapter_for_pe \
  --disable_quality_filtering \
  --length_required 30 \
  --thread 8 \
  --html qc_reports_trimmed/SRR11998467_fastp.html \
  --json qc_reports_trimmed/SRR11998467_fastp.json 

#Starvation
fastp \
  -i starvation/SRR11998473_1.fastq \
  -I starvation/SRR11998473_2.fastq \
  -o trimmed/SRR11998473_1.clean.fastq.gz \
  -O trimmed/SRR11998473_2.clean.fastq.gz \
  --detect_adapter_for_pe \
  --disable_quality_filtering \
  --length_required 30 \
  --thread 8 \
  --html qc_reports_trimmed/SRR11998473_fastp.html \
  --json qc_reports_trimmed/SRR11998473_fastp.json 


```

```bash
cd ~/groupproject
mkdir qc_reports/trimmed

fastqc -t 8 trimmed/*.fastq.gz -o qc_reports_trimmed

cd qc_reports_trimmed
multiqc . -o multiqc_trimmed

```
- **Tasks Performed:** 
  - Ran `fastp` on all paired-end reads for each condition
  - Removed adapter sequences from raw reads. Quality filtering was disabled due to masked quality scores.
  - Generated per-sample trimming reports (`.html` and `.json`)    
  - Reran `fastqc` on the trimmed files
  - Ran `multiqc` on the trimmed files

- **Deliverables:**
### 1-2 QC plots for cleaned data
  - [MultiQC Report]()
  - [Adapter Content Plot]()
  - [Sequence Duplication Levels]()
  - [Per Sequence GC Content]()

### Cleaned Read Counts and Duplication Rates (from multiqc_general_stats.txt)
| Sample        | fastp Dup % | fastp Q30 Rate | fastp Q30 Bases | fastp GC % | % Surviving | % Adapters | FastQC Dup % | FastQC GC % | Avg Length (bp) | Total Reads |
|---------------|-------------|----------------|-----------------|------------|-------------|------------|--------------|-------------|-----------------|-------------|
| SRR11998457_1 | 66.3%       | 0.0            | 0.0             | 50.5%      | 92.6%       | 7.9%       | 88.4%        | 50.0%       | 141.2           | 21.0 M      |
| SRR11998457_2 | —           | —              | —               | —          | —           | —          | 89.9%        | 50.0%       | 149.9           | 21.0 M      |
| SRR11998467_1 | 66.7%       | 0.0            | 0.0             | 51.2%      | 92.0%       | 6.3%       | 86.4%        | 51.0%       | 141.8           | 9.4 M       |
| SRR11998467_2 | —           | —              | —               | —          | —           | —          | 87.5%        | 50.0%       | 150.7           | 9.4 M       |
| SRR11998473_1 | 67.8%       | 0.0            | 0.0             | 50.9%      | 90.6%       | 9.8%       | 86.9%        | 51.0%       | 141.2           | 4.8 M       |
| SRR11998473_2 | —           | —              | —               | —          | —           | —          | 87.4%        | 50.0%       | 149.8           | 4.8 M       |

### Raw vs. cleaned read counts.
| Sample        | Condition   | Raw Reads | Cleaned Reads | % Surviving |
|---------------|-------------|-----------|---------------|-------------|
| SRR11998457   | Acidic      | 22.7 M    | 21.0 M        | 92.6%       |
| SRR11998467   | Oxidative   | 10.2 M    | 9.4 M         | 92.0%       |
| SRR11998473   | Starvation  | 5.4 M     | 4.8 M         | 90.6%       |



- **Interpretation:**  
  - Read counts: Read count less than raw reads (22.7M -> 21.0M, 10.2M -> 9.4M, 5.4M -> 4.8M). ~90% of sequences were passed through cleaning. Low quality bases were not trimmed due to hidden scores. This drop is due to adapter removal and fixed trimming.
  - Duplication: 66-68% in fastp, 86-90% in fastqc. Still remains high but consistent with RNA-seq.
  - GC content: Still stable at 50-51%. Trimming did not change base composition. 
  - Read length: Varied between 141-151bp due to adapter trimming clipping reads to different lengths.

## Download Salmonella Reference Genome

- **Tasks Performed:** 
  - Confirmed SRR11998457, SRR11998467, and SRR11998473 in BioProject PRJNA638918 are from Salmonella enterica Typhimurium
  - Found reference genome for Salmonella enterica serovar Typhimurium LT2 (GCF_000006945.2 (ASM694v2))
  - Downloaded genome FASTA (.fna) and annotation file (.gff)
  - Renamed files to styphimurium.fna and styphimurium.gff and moved into a reference directory

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz
gunzip *.gz

mv GCF_000006945.2_ASM694v2_genomic.fna styphimurium.fna
mv GCF_000006945.2_ASM694v2_genomic.gff styphimurium.gff

mkdir reference

mv styphimurium.* reference/

```
## Alignment
- Tools Used:
  - HISAT2
    

- **Tasks Performed:** 
  - 
  - 

- **Deliverables:**
  - Plots of alignment metrics (total reads, mapping %)
  - 

## Annotation and Quantification
- Tools Used: Salmon
  - featureCounts
    
- **Tasks Performed:** 
  - Setup + Folders
  - #install Salmon via conda
mamba create -n salmon-env -c bioconda -c conda-forge salmon pigz parallel -y
conda activate salmon-env

  - 

- **Deliverables:**
  - Top 10 expressed genes table
  - Top 10 expressed heatmap
  - Functional enrichment
