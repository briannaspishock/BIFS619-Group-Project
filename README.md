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
  - Subread (https://subread.sourceforge.net/) for feature counting
```bash
#install packages
sudo apt install samtools
sudo apt install hisat2
sudo apt install subread

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


```

- **Tasks Performed:** 
  - 
  - 

- **Deliverables:**
- 1-2 tables with alignment metrics
  - Identify top 10 genes expressed
  - Optional: perform functional enrichment 

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
mkdir -p quant tables logs
```

  -Get more reference for Salmonella

```bash
cd reference
#Genome coding sequences (CDS), proteins sequence
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_cds_from_genomic.fna.gz

#decompress all the .gz files
pigz -d *.gz

#Renamed files
mv GCF_000006945.2_ASM694v2_cds_from_genomic.fna styphimurium_cds_from_genomic.fna

cd ..

```

  -Build Salmon index (from CDS)
```bash
salmon index \
  -t reference/styphimurium_cds_from_genomic.fna \
-i reference/salmon_cds_index
```

  -Quantify each condition with Salmon
```bash
#take the raw RNA-seq reads (*_1.fastq and *_2.fastq) from each folder and produces expression quantification results (TPM, counts, etc.) for each condition.
for cond in acidic oxidative starvation; do
  r1=$(ls ${cond}/*_1.fastq)
  r2=$(ls ${cond}/*_2.fastq) 
  salmon quant \   
    -i reference/salmon_cds_index \
    -l A \
    -1 "$r1" -2 "$r2" \
    -p 8 --gcBias --validateMappings \
    -o quant/${cond} \
    2> logs/${cond}.salmon.log
Done
```

  -Check of mapping rates
```bash
grep -E "Mapping rate|chosen|libType" logs/*.salmon.log
```

**Annotation:**

  -Map transcript/CDS IDs → locus_tag from the CDS FASTA headers
```bash
grep "^>" reference/styphimurium_cds_from_genomic.fna \
| sed 's/^>//' \
| awk '{
  tx=$1;
  lt="";
  # accept both bracketed and unbracketed forms
  if (match($0,/locus_tag=([A-Za-z0-9_\.-]+)/,m)) lt=m[1];
  else if (match($0,/\[locus_tag=([A-Za-z0-9_\.-]+)\]/,m)) lt=m[1];
  if (lt!="") print tx"\t"lt;
}' > tables/tx2gene.tsv
```

  -build a locus_tag → gene_name/product table from the GFF
```bash
awk -F'\t' '$3=="gene"{
  split($9,a,";");
  lt=""; gn=""; pr="";
  for(i in a){
    if(a[i] ~ /^locus_tag=/) lt=substr(a[i],11);
    else if(a[i] ~ /^gene=/)     gn=substr(a[i],6);
    else if(a[i] ~ /^Name=/)     pr=substr(a[i],6);
    else if(a[i] ~ /^product=/)  pr=substr(a[i],9);
  }
  if(lt!="") print lt"\t" (gn==""?"-":gn) "\t" (pr==""?"-":pr);
}' reference/styphimurium.gff \
> tables/gene_annot.tsv
```
  -Top-10 expressed genes by TPM per condition
```bash
# Sort tx2gene.tsv by the first column (transcript ID)
sort -k1,1 tables/tx2gene.tsv -o tables/tx2gene.tsv

for cond in acidic oxidative starvation; do
  # Extract transcript and TPM
  awk 'NR>1{print $1"\t"$4}' quant/${cond}/quant.sf \
    | sort -k1,1 > tables/${cond}.tx_tpm.tsv

  # Sum TPM per gene
  join -t $'\t' -1 1 -2 1 tables/${cond}.tx_tpm.tsv tables/tx2gene.tsv \
    | awk -F'\t' '{tpm[$3]+=$2} END{for(g in tpm) print g"\t"tpm[g]}' \
    | sort -k2,2nr | head -n 10 > tables/top10_${cond}_TPM.raw.tsv

  # Annotate gene names and products
  sort -k1,1 tables/top10_${cond}_TPM.raw.tsv > tables/a.tmp
  sort -k1,1 tables/gene_annot.tsv > tables/b.tmp
  join -t $'\t' -1 1 -2 1 tables/a.tmp tables/b.tmp \
    | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,($3==""?"-":$3),($4==""?"-":$4)}' \
    > tables/top10_${cond}_TPM.annot.tsv

sort -k2,2nr tables/top10_${cond}_TPM.annot.tsv -o tables/top10_${cond}_TPM.annot.tsv

  rm -f tables/a.tmp tables/b.tmp
  echo "Rebuilt clean top 10 for $cond -> tables/top10_${cond}_TPM.annot.tsv"
done


#Add headers to all 3 files
for cond in acidic oxidative starvation; do
  file="tables/top10_${cond}_TPM.annot.tsv"
  tmp="tables/top10_${cond}_TPM.annot.tmp"
  echo -e "locus_tag\tTPM_sum\tgene_name\tproduct" > "$tmp"
  cat "$file" >> "$tmp"
  mv "$tmp" "$file"
  echo "Added header to $file"
done
```

  -Each
  
  tables/top10_acidic_TPM.annot.tsv

   tables/top10_oxidative_TPM.annot.tsv
   
   tables/top10_starvation_TPM.annot.tsv
   
  contains:
    
  locus_tag : The unique gene ID in the genome (from the annotation .gff file)
  
  TPM_sum : The total expression level for that gene, summed across its transcripts. The higher the TPM, the more abundant the RNA
  
  gene_name : The standard gene symbol
  
  product : The protein or functional product encoded by the gene
    
   

- **Deliverables:**
  - Top 10 expressed genes table
  - Top 10 expressed heatmap
  - Functional enrichment
