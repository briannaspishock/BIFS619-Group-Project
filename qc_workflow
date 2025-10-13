#!/usr/bin/env bash
# qc_workflow.sh
# BIFS619 Group Project Salmonella QC workflow, Brianna Spishock
# - download Salmonella RNA-seq data (3 runs) into ~/groupproject/{acidic,oxidative,starvation}
# - QC (FastQC + MultiQC) on raw reads
# - cleaning (fastp) + QC (FastQC + MultiQC) on trimmed reads
# - download Salmonella reference genome + annotation into ~/groupproject/reference

set -euo pipefail

echo "=== downloading Salmonella RNA-seq data ==="

# enter analysis directory
cd ~/groupproject

# install sratool-toolkit to use prefetch and other tools
sudo apt update
sudo apt install -y sra-toolkit

# acidic stress 
mkdir -p acidic
prefetch SRR11998457
fasterq-dump SRR11998457 --split-files --threads 8 --outdir acidic

# oxidative stress
mkdir -p oxidative
prefetch SRR11998467
fasterq-dump SRR11998467 --split-files --threads 8 --outdir oxidative

# starvation stress
mkdir -p starvation
prefetch SRR11998473
fasterq-dump SRR11998473 --split-files --threads 8 --outdir starvation

echo "=== STEP 2: QC with fastqc, MultiQC ==="

# install fastqc + multiqc (fastqc is needed before next step)
sudo apt install -y fastqc
# we'll install multiqc now or right before use; either is fine
sudo apt install -y multiqc

# create an output folder for FastQC reports
mkdir -p qc_reports

# run FastQC on all FASTQ files inside each condition folder
fastqc acidic/*.fastq -o qc_reports
fastqc oxidative/*.fastq -o qc_reports
fastqc starvation/*.fastq -o qc_reports

# put all reports together using MultiQC
cd qc_reports
multiqc .
# generate QC plots 
multiqc . --export
cd ..

echo "=== STEP 3: read cleaning with fastp, MultiQC ==="

# install fastp
sudo apt install -y fastp

# make directory for reports
cd ~/groupproject
mkdir -p trimmed
mkdir -p qc_reports_trimmed

# acidic
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

# oxidative
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

# starvation
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

# fastQC on trimmed files
fastqc -t 8 trimmed/*.fastq.gz -o qc_reports_trimmed

# multiQC on trimmed files
cd qc_reports_trimmed
multiqc . -o multiqc_trimmed
cd ..

echo "=== STEP 4: downloading Salmonella reference genome ==="

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz
gunzip *.gz

mv GCF_000006945.2_ASM694v2_genomic.fna styphimurium.fna
mv GCF_000006945.2_ASM694v2_genomic.gff styphimurium.gff

mkdir -p reference
mv styphimurium.* reference/

echo "------------------------------------------------------------"
echo "QC completed!"
echo "raw reads:        ~/groupproject/{acidic,oxidative,starvation}"
echo "raw QC:           ~/groupproject/qc_reports (includes MultiQC + exports)"
echo "trimmed reads:    ~/groupproject/trimmed"
echo "trimmed QC:       ~/groupproject/qc_reports_trimmed (multiqc_trimmed/)"
echo "reference:        ~/groupproject/reference/{styphimurium.fna, styphimurium.gff}"
echo "------------------------------------------------------------"
