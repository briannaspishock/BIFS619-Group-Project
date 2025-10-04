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
  - QC plots
- [MultiQC Report](multiqc_report_1.html)
- [Adapter Content Plot](multiqc_plots/mqc_fastqc_adapter_content_plot_1.png)
- [Sequence Duplication Levels](multiqc_plots/mqc_fastqc_sequence_duplication_levels_plot_1.png)
- [Per Sequence GC Content](multiqc_plots/mqc_fastqc_per_sequence_gc_content_plot_Percentages.png)

Insert Table
### Raw Read Counts and Duplication Rates

| Sample        | Condition   | Raw Reads | Duplication % | GC % | Avg Length (bp) |
|---------------|-------------|-----------|---------------|------|-----------------|
| SRR11998457_1 | Acidic      | 22.7 M    | 100%          | 50%  | 143             |
| SRR11998457_2 | Acidic      | 22.7 M    | 100%          | 50%  | 151             |
| SRR11998467_1 | Oxidative   | 10.2 M    | 100%          | 51%  | 143             |
| SRR11998467_2 | Oxidative   | 10.2 M    | 100%          | 51%  | 151             |
| SRR11998473_1 | Starvation  | 5.4 M     | 100%          | 51%  | 143             |
| SRR11998473_2 | Starvation  | 5.4 M     | 100%          | 50%  | 151             |

 

- **Interpretation:**  
- Read counts: Range from 22.7M in acidic to 5.4M in starvation. This reflects sequencing depth across different conditions.
- Duplication: Reported as 100% in all samples.
- GC content: Stable at 50-51%
- Read length: forward reads are 143bp and reverse reads are 151bp, consistent with expected paired-end design.

## Read Cleaning
- Tools Used:
- Fastp
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
  - Removed low-quality bases and adapter sequences from raw reads  
  - Generated per-sample trimming reports (`fastp_report.html` and `fastp_report.json`)    
  - Reran 'fastqc' on the trimmed files

- **Deliverables:**
  - QC plots for cleaned data


Insert Table
### Cleaned Read Counts and Duplication Rates

| Sample               | Condition   | Clean Reads | Duplication % | GC % | Avg Length (bp) |
|----------------------|-------------|-------------|---------------|------|-----------------|
| SRR11998457_1.clean  | Acidic      | 21.0 M      | 100%          | 50%  | 30–143          |
| SRR11998457_2.clean  | Acidic      | 21.0 M      | 100%          | 50%  | 31–151          |
| SRR11998467_1.clean  | Oxidative   | 9.4 M       | 100%          | 51%  | 30–143          |
| SRR11998467_2.clean  | Oxidative   | 9.4 M       | 100%          | 50%  | 31–151          |
| SRR11998473_1.clean  | Starvation  | 4.8 M       | 100%          | 51%  | 31–143          |
| SRR11998473_2.clean  | Starvation  | 4.8 M       | 100%          | 50%  | 31–151          |

- **Interpretation:**  
- Read counts: Read count less than raw reads (22.7M -> 21.0M, 10.2M -> 9.4M, 5.4M -> 4.8M). Approximately 90% of sequences were passed through cleaning. Low quality bases were no trimmed due to hidden scores. This drop is due to adapter removal and fixed trimming.
- Duplication: Still showing 100%. 
- GC content: Still stable at 50-51%. Trimming did not change base composition. 
- Read length: Varied between 30-151bp due to adapter trimming clipping reads to different lengths.

