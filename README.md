# RNA-Seq Analysis Project
This repository contains code, results, and documentation for our BIFS619 group project. The goal of this project is to analyze microbial RNA-seq data using a reproducible pipeline, from raw FASTQ files to quality control, alignment, quantification, and visualization.

/data       → placeholder for raw and processed data (not uploaded; SRR IDs only)  
/scripts    → all code (QC, trimming, alignment, quantification, plotting)  
/results    → figures, tables, plots generated from analysis  
/reports    → Analysis Report, Team Report, Individual Reports  

# Workflow Overview
Data Retrieval – Download FASTQ files from SRA (3 samples).
Quality Control – FastQC + MultiQC reports for raw reads.
Read Cleaning – Adapter/quality trimming using fastp or Cutadapt/Trimmomatic.
Alignment – Reads mapped to reference genome (HISAT2/Salmon).
Quantification – Expression levels quantified (Salmon/featureCounts).
Analysis & Visualization – Heatmaps, top expressed genes, and biological interpretation.
Reporting – Analysis Report + Team Report compiled with results and discussion.

# Tools & Dependencies
FastQC / MultiQC
fastp 
HISAT2 
featureCounts
Python (matplotlib, seaborn) for plots
Environment: managed with Conda for reproducibility

# Deliverables
Quality control plots + summary table
Pre/post trimming comparison
Alignment metrics table/plots
Expression table (TPM/CPM)
Heatmap of top 10–20 expressed genes
Final reports (Analysis + Team + Individual)
