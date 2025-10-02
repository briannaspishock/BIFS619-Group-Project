# Group Project: RNA-Seq Workflow (Microbial)

## Preparation
- **GitHub Setup**
  - Team members: Brianna Spishock, Youngju Jeon, Brendan O'Brien 

## Data Retrieval
- **NCBI SRA accessions used:**
  - SRR11998473 -starvation
  - SRR11998467 -oxidative
  - SRR11998457 -acidic

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
```bash
# Create an output folder for FastQC reports
mkdir -p qc_reports

# Run FastQC on all FASTQ files inside each condition folder
fastqc acidic/*.fastq -o qc_reports
fastqc oxidative/*.fastq -o qc_reports
fastqc starvation/*.fastq -o qc_reports
```
  - MultiQC

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
Insert QC plots
- [MultiQC Report](multiqc_report_1.html)
- [Adapter Content Plot](multiqc_plots/mqc_fastqc_adapter_content_plot_1.png)
- [Sequence Duplication Levels](multiqc_plots/mqc_fastqc_sequence_duplication_levels_plot_1.png)
- [Per Sequence GC Content](multiqc_plots/mqc_fastqc_per_sequence_gc_content_plot_Percentages.png)

Insert Table
### Raw Read Counts and Duplication Rates

| Sample       | Condition  | Raw Reads   | Duplication % | GC % | Avg Length (bp) |
|--------------|------------|-------------|---------------|------|-----------------|
| SRR11998457_1 | Acidic     | 4.5 M       | 76%           | 50%  | 143 |
| SRR11998457_2 | Acidic     | 4.5 M       | 79%           | 50%  | 151 |
| SRR11998467_1 | Oxidative  | 10.2 M      | 88%           | 51%  | 143 |
| SRR11998467_2 | Oxidative  | 10.2 M      | 88%           | 51%  | 151 |
| SRR11998473_1 | Starvation | 5.4 M       | 88%           | 51%  | 143 |
| SRR11998473_2 | Starvation | 5.4 M       | 88%           | 50%  | 151 |
 

- **Interpretation:**  
- Raw read counts were highest in the oxidative condition (~20M), moderate in starvation (~10M), and lowest in acidic (~9M). GC content was stable at ~50–51%. High duplication levels (>75%) in RNA-seq raw reads. There is a need for trimming and filtering in the read cleaning step.

## Read Cleaning
- Tools Used:
- Fastp
 ```bash
#Acidic
fastp -i acidic/SRR11998457_1.fastq \
      -I acidic/SRR11998457_2.fastq \
      -o acidic/SRR11998457_1.trimmed.fastq \
      -O acidic/SRR11998457_2.trimmed.fastq \
      -h acidic/fastp_report.html \
      -j acidic/fastp_report.json

#Oxidative
fastp \
  -i oxidative/SRR11998467_1.fastq \
  -I oxidative/SRR11998467_2.fastq \
  -o oxidative/SRR11998467_1.trimmed.fastq \
  -O oxidative/SRR11998467_2.trimmed.fastq \
  -h oxidative/fastp_report.html \
  -j oxidative/fastp_report.json

#Starvation
fastp \
  -i starvation/SRR11998473_1.fastq \
  -I starvation/SRR11998473_2.fastq \
  -o starvation/SRR11998473_1.trimmed.fastq \
  -O starvation/SRR11998473_2.trimmed.fastq \
  -h starvation/fastp_report.html \
  -j starvation/fastp_report.json


```
- **Tasks Performed:** 
  - Ran `fastp` on all paired-end reads for each condition
  - Removed low-quality bases and adapter sequences from raw reads  
  - Generated per-sample trimming reports (`fastp_report.html` and `fastp_report.json`)    
  - Reran 'fastqc' on the trimmed files
  - 
```bash
cd ~/groupproject
mkdir qc_reports/trimmed

# Acidic
fastqc acidic/*trimmed.fastq -o qc_reports/trimmed

# Oxidative
fastqc oxidative/*trimmed.fastq -o qc_reports/trimmed

# Starvation
fastqc starvation/*trimmed.fastq -o qc_reports/trimmed
# Acidic
fastqc acidic/*trimmed.fastq -o qc_reports/trimmed

# Oxidative
fastqc oxidative/*trimmed.fastq -o qc_reports/trimmed

# Starvation
fastqc starvation/*trimmed.fastq -o qc_reports/trimmed

# Rerun MultiQC on trimmed files
cd qc_reports/trimmed
multiqc . -o ../multiqc_trimmed
```
- **Deliverables:**  
- **Interpretation:**  
