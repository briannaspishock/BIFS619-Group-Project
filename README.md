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
  -Setup + Folders

```bash
#install Salmon via conda
mamba create -n salmon-env -c bioconda -c conda-forge salmon pigz parallel -y
conda activate salmon-env
cd ~/groupproject
mkdir -p ref quant tables logs
```

  -Get reference (CDS + GFF) for Salmonella

```bash
cd ref
#Genome coding sequences (CDS), proteins, and GFF annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_protein.faa.gz

#decompress all the .gz files
pigz -d *.gz
cd ..

```

  -Build Salmon index (from CDS)
```bash
salmon index \
  -t ref/GCF_000006945.2_ASM694v2_cds_from_genomic.fna \
-i ref/salmon_cds_index
```

  -Quantify each condition with Salmon
```bash
#take the raw RNA-seq reads (*_1.fastq and *_2.fastq) from each folder and produces expression quantification results (TPM, counts, etc.) for each condition.
for cond in acidic oxidative starvation; do
#find the forward read file (_1.fastq) inside the folder for that condition and stores the file name in a variable called r1
  r1=$(ls ${cond}/*_1.fastq)
#does the same for the reverse read file (_2.fastq)
  r2=$(ls ${cond}/*_2.fastq) 
#calls Salmon, a tool that estimates how much each gene (or transcript) is expressed by aligning reads to the reference CDS database (from salmon_cds_index)
  salmon quant \   
#Tells Salmon which index (the pre-built reference from your CDS FASTA) to use
    -i ref/salmon_cds_index \
#Tells Salmon to automatically detect the library type
    -l A \
#Specifies the paired-end read files for quantification
    -1 "$r1" -2 "$r2" \
#Use 8 CPU threads to speed up processing. Enable GC-bias correction, which adjusts expression estimates if GC-rich genes are over- or under-represented in sequencing.
    -p 8 --gcBias --validateMappings \
#Writes all the results into a separate output folder
    -o quant/${cond} \
#Redirects all log messages and errors from Salmon into a log file per condition
    2> logs/${cond}.salmon.log
Done
```

  -Check of mapping rates
```bash
grep -E "Mapping rate|chosen|libType" logs/*.salmon.log
```

  -Map transcript/CDS IDs → locus_tag from the CDS FASTA headers
```bash
grep "^>" ref/GCF_000006945.2_ASM694v2_cds_from_genomic.fna \
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
}' ref/GCF_000006945.2_ASM694v2_genomic.gff \
> tables/gene_annot.tsv
```
  -Top-10 expressed genes by TPM per condition
```bash
# Sort tx2gene.tsv by the first column (transcript ID)
sort -k1,1 tables/tx2gene.tsv -o tables/tx2gene.tsv

for cond in acidic oxidative starvation; do
  awk 'NR>1{print $1"\t"$4}' quant/${cond}/quant.sf \
    | sort -k1,1 > tables/${cond}.tx_tpm.tsv

  join -t $'\t' -1 1 -2 1 tables/${cond}.tx_tpm.tsv tables/tx2gene.tsv \
    | awk -F'\t' '{tpm[$3]+=$2} END{for(g in tpm) print g"\t"tpm[g]}' \
    | sort -k2,2nr | head -n 10 > tables/top10_${cond}_TPM.raw.tsv

  sort -k1,1 tables/top10_${cond}_TPM.raw.tsv > tables/a.tmp
  sort -k1,1 tables/gene_annot.tsv > tables/b.tmp
  join -t $'\t' -1 1 -2 1 tables/a.tmp tables/b.tmp \
    | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,($3==""?"-":$3),($4==""?"-":$4)}' \
    > tables/top10_${cond}_TPM.annot.tsv
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
  echo "✅ Added header to $file"
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
