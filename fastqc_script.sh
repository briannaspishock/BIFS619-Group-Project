#!/usr/bin/env bash
# quick FASTQC on downloaded SRA files
# bifs619 group 5

set -euo pipefail
OUTDIR=~/sra_work/fastqc_reports
mkdir -p "$OUTDIR"


#run fastqc on gzipped files
fastqc -t 8 ~/sra_work/fastq/*.fastq.gz -o "$OUTDIR"
