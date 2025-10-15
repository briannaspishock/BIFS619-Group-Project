# Group Project: RNA-Seq Workflow (Microbial)

## Preparation
- **GitHub Setup**
  - Team members: Brianna Spishock (QC), Youngju Jeon (Annotation), Brendan O'Brien (Alignment)

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

# install sratool-toolkit to use prefetch and other tools
sudo apt install sra-toolkit

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
#install multiqc
sudo apt install multiqc
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
  - QC plots in analysis file

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
  - MultiQC

 ```bash
# install fastp
sudo apt install fastp

# make directory for reports
cd ~/groupproject
mkdir trimmed
mkdir qc_reports_trimmed

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
  - QC plots in analysis doc

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
  - Downloaded genome FASTA (.fna), CDS FASTA(cds_from_genomic.fna) and annotation file (.gff)
  - Renamed files to styphimurium.fna, styphimurium.gff and styphimurium_cds_from_genomic.fna then moved into a reference directory

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_cds_from_genomic.fna.gz
gunzip *.gz

mv GCF_000006945.2_ASM694v2_genomic.fna styphimurium.fna
mv GCF_000006945.2_ASM694v2_genomic.gff styphimurium.gff
mv GCF_000006945.2_ASM694v2_cds_from_genomic.fna styphimurium_cds_from_genomic.fna

mkdir reference

mv styphimurium.* reference/

```
### The full quality control workflow is provided in [`qc_workflow.sh`](./qc_workflow).

This script:

- Downloads Salmonella RNA-seq data (acidic, oxidative, starvation conditions)
- Runs `FastQC` and `MultiQC` on raw reads
- Cleans reads with `fastp` and re-runs QC on trimmed reads
- Downloads the Salmonella reference genome and annotation files

### Usage

```bash
#make executable
chmod +x qc_workflow.sh

#run workflow
./qc_workflow.sh
```

## Alignment
- Tools Used:
  - HISAT2
  - Samtools
  - Samtools flagstat for alignment data
```bash
#install packages
sudo apt install samtools
sudo apt install hisat2

#change directory to ref genome
cd reference/

#create index files
hisat2-build styphimurium.fna styph_index

# go back to project folder
cd ../

#create alignment directory
mkdir -p alignment_bam

# perform alignments with Hisat2  using cleaned fastq (these may take a few minutes)
#acidic
hisat2 -p 8 -x reference/styph_index \
 -1 trimmed/SRR11998457_1.clean.fastq.gz \
 -2 trimmed/SRR11998457_2.clean.fastq.gz | \
 samtools view -b - | \
 samtools sort -o alignment_bam/acidic_sorted.bam

#oxidative
hisat2 -p 8 -x reference/styph_index \
 -1 trimmed/SRR11998467_1.clean.fastq.gz \
 -2 trimmed/SRR11998467_2.clean.fastq.gz | \
 samtools view -b - | \
 samtools sort -o alignment_bam/oxidative_sorted.bam

#starvation
hisat2 -p 8 -x reference/styph_index \
 -1 trimmed/SRR11998473_1.clean.fastq.gz \
 -2 trimmed/SRR11998473_2.clean.fastq.gz | \
 samtools view -b - | \
 samtools sort -o alignment_bam/starvation_sorted.bam


# index the acidic sample
samtools index alignment_bam/acidic_sorted.bam

# Index the oxidative sample
samtools index alignment_bam/oxidative_sorted.bam

# Index the starvation sample
samtools index alignment_bam/starvation_sorted.bam

#alignment statistics
 mkdir -p alignment_stats

# Acidic sample
samtools flagstat alignment_bam/acidic_sorted.bam > alignment_stats/acidic_stats.txt

# Oxidative sample
samtools flagstat alignment_bam/oxidative_sorted.bam > alignment_stats/oxidative_stats.txt

# Starvation sample
samtools flagstat alignment_bam/starvation_sorted.bam > alignment_stats/starvation_stats.txt

```

- **Tasks Performed:** 
  - created index files from reference genome
  - created sam files from trimmed, cleaned fastq 
  - align, converted sam to bam
  - index bam files
  - used flagstat for alignment metrics  

- **Deliverables:**
- 1-2 tables with alignment metrics

  | Sample Name | Condition   | Total Reads  | Mapping % |
| :---------- | :---------- | :----------- | :-------- |
| SRR11998457 | Acidic      |  43363202    |  95.72%     |
| SRR11998467 | Oxidative   |  19077236    |  95.94% |
| SRR11998473 | Starvation  |  9994445     | 95.74% |

## Annotation and Quantification
- Tools Used: Salmon
  - featureCounts
    
- **Tasks Performed:** 

**Quantification:**

  -Setup + Folders

```bash
#install Salmon via conda
mamba create -n salmon-env -c bioconda -c conda-forge salmon pigz parallel -y
conda activate salmon-env
cd ~/groupproject
mkdir -p quant tables maps
```


  -Build Salmon index (from CDS)
```bash
salmon index \
  -t reference/styph_cds.fna \
  -i reference/salmon_cds_index \
  --kmerLen 31
```

  -Quantify each condition with Salmon
```bash
#take the raw RNA-seq reads (*_1.fastq and *_2.fastq) from each folder and produces expression quantification results (TPM, counts, etc.) for each condition.
# Acidic
salmon quant -i reference/salmon_cds_index -l A \
  -1 trimmed/SRR11998457_1.clean.fastq.gz \
  -2 trimmed/SRR11998457_2.clean.fastq.gz \
  -p 8 --validateMappings --gcBias --seqBias \
  -o quant/acidic

# Oxidative
salmon quant -i reference/salmon_cds_index -l A \
  -1 trimmed/SRR11998467_1.clean.fastq.gz \
  -2 trimmed/SRR11998467_2.clean.fastq.gz \
  -p 8 --validateMappings --gcBias --seqBias \
  -o quant/oxidative

# Starvation
salmon quant -i reference/salmon_cds_index -l A \
  -1 trimmed/SRR11998473_1.clean.fastq.gz \
  -2 trimmed/SRR11998473_2.clean.fastq.gz \
  -p 8 --validateMappings --gcBias --seqBias \
  -o quant/starvation
```

  -Build a robust annotation map (CDS → protein_id → gene/product)
```bash
# FASTA→protein map
awk '
  /^>/ {
    tid=$1; sub(/^>/,"",tid);
    pid="";
    if (match($0, /protein_id=([^]]+)/, m)) { 
      pid=m[1];
    }
    if (pid!="") print tid "\t" pid;
  }
' reference/styphimurium_cds_from_genomic.fna > maps/fastaid2protein.tsv

# GFF protein→gene/product map
awk -F'\t' '
$3=="CDS" {
  pid=""; ltag=""; prod="";
  split($9, a, ";");
  for (i in a) {
    split(a[i], kv, "=");
    if (kv[1]=="protein_id") pid=kv[2];
    else if (kv[1]=="locus_tag") ltag=kv[2];
    else if (kv[1]=="product")   prod=kv[2];
  }
  if (pid!="") print pid"\t"ltag"\t"prod;
}' reference/styphimurium.gff > maps/protein2gene_product.tsv

#Sort and join on protein_id
sort -k2,2 maps/fastaid2protein.tsv > maps/fastaid2protein.sorted.tsv   # key is col2
sort -k1,1 maps/protein2gene_product.tsv > maps/protein2gene_product.sorted.tsv  # key is col1

# Join on protein_id, then reorder to (transcript_id, gene, product)
join -t $'\t' -1 2 -2 1 maps/fastaid2protein.sorted.tsv maps/protein2gene_product.sorted.tsv \
| awk -v OFS='\t' '{tid=$2; gene=$3; prod=$4; print tid, gene, prod}' \
> maps/tx2gene_product.tsv
```

  -Top-10 tables
```bash
#Gene-level tables from Salmon quant.sf
#Acidic (Top 10 TPM)
awk -F'\t' 'NR==FNR{m[$1]=$2"\t"$3; next} FNR>1{print $4"\t"m[$1]"\t"$1"\t"$5}' \
  maps/tx2gene_product.tsv \
  quant/acidic/quant.sf \
| sort -nrk1 \
| head -10 \
| awk 'BEGIN{OFS="\t"; print "TPM","Gene","Product","CDS_ID","NumReads"} {print}' \
> tables/top10_acidic.tsv

#Oxidative (Top 10 TPM)
awk -F'\t' 'NR==FNR{m[$1]=$2"\t"$3; next} FNR>1{print $4"\t"m[$1]"\t"$1"\t"$5}' \
  maps/tx2gene_product.tsv \
  quant/oxidative/quant.sf \
| sort -nrk1 | head -10 \
| awk 'BEGIN{OFS="\t"; print "TPM","Gene","Product","CDS_ID","NumReads"} {print}' \
> tables/top10_oxidative.tsv

#Starvation (Top 10 TPM)
awk -F'\t' 'NR==FNR{m[$1]=$2"\t"$3; next} FNR>1{print $4"\t"m[$1]"\t"$1"\t"$5}' \
  maps/tx2gene_product.tsv \
  quant/starvation/quant.sf \
| sort -nrk1 | head -10 \
| awk 'BEGIN{OFS="\t"; print "TPM","Gene","Product","CDS_ID","NumReads"} {print}' \
> tables/top10_starvation.tsv

```

  -Each
  
tables/top10_acidic.tsv

tables/top10_oxidative.tsv

tables/top10_starvation.tsv

 

- **Deliverables:**
  - Top 10 expressed genes table
  - Top 10 expressed heatmap
  - Functional enrichment
