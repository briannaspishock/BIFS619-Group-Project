#!/usr/bin/env bash
# trimming + QC for paired-end FASTQs
# bifs619 group project group 5

set -euo pipefail

# ================== CONFIG ==================
RAW_DIR="${RAW_DIR:-$HOME/sra_work/fastq}"            # where raw FASTQs are
TRIM_DIR="${TRIM_DIR:-$HOME/sra_work/trimmed}"        # where trimmed outputs go
QC_DIR="${QC_DIR:-$HOME/sra_work/fastqc_reports_trimmed}"  # FastQC/MultiQC outputs
THREADS="${THREADS:-8}"                                # adjust to your VM cores
# ============================================

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 SRR123 SRR456 ..."
  exit 1
fi

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$QC_DIR"

need() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: '$1' not found in PATH." >&2
    exit 1
  fi
}
need fastp
need fastqc

# helper: choose .fastq.gz if present, else .fastq
pick_fastq() {
  local base="$1"
  if [[ -f "${base}.fastq.gz" ]]; then
    printf "%s.fastq.gz" "$base"
  elif [[ -f "${base}.fastq" ]]; then
    printf "%s.fastq" "$base"
  else
    printf ""
  fi
}

# ================== TRIM PER SRR ==================
for srr in "$@"; do
  echo "=== fastp: $srr ==="
  r1_base="$RAW_DIR/${srr}_1"
  r2_base="$RAW_DIR/${srr}_2"

  R1="$(pick_fastq "$r1_base")"
  R2="$(pick_fastq "$r2_base")"
  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: Could not find inputs for $srr" >&2
    echo "Looked for: ${r1_base}.fastq.gz|.fastq and ${r2_base}.fastq.gz|.fastq" >&2
    exit 1
  fi

  out1="$TRIM_DIR/${srr}_1.trimmed.fastq.gz"
  out2="$TRIM_DIR/${srr}_2.trimmed.fastq.gz"
  html="$TRIM_DIR/${srr}.fastp.html"
  json="$TRIM_DIR/${srr}.fastp.json"

  fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$out1" \
    -O "$out2" \
    -w "$THREADS" \
    -h "$html" \
    -j "$json"

  echo "✓ Trimmed $srr → $(basename "$out1"), $(basename "$out2")"
done

# ================== FASTQC ON TRIMMED ==================
echo "=== FastQC on trimmed reads ==="
fastqc -t "$THREADS" "$TRIM_DIR"/*trimmed.fastq.gz -o "$QC_DIR"

# ================== MULTIQC (optional) ==================
if command -v multiqc >/dev/null 2>&1; then
  echo "=== MultiQC summary ==="
  (cd "$QC_DIR" && multiqc . -o .)
  echo "✓ MultiQC report: $QC_DIR/multiqc_report.html"
else
  echo "NOTE: 'multiqc' not found; skipping MultiQC. Install with 'conda install -c bioconda multiqc' if you want a combined report."
fi

# ================== SIMPLE READ COUNT SUMMARY ==================
echo -e "SRR\tR1_reads(trimmed)\tR2_reads(trimmed)" > "$TRIM_DIR/trimmed_read_counts.tsv"
for srr in "$@"; do
  r1_lines=$(zcat "$TRIM_DIR/${srr}_1.trimmed.fastq.gz" | wc -l)
  r2_lines=$(zcat "$TRIM_DIR/${srr}_2.trimmed.fastq.gz" | wc -l)
  echo -e "${srr}\t$((r1_lines/4))\t$((r2_lines/4))" >> "$TRIM_DIR/trimmed_read_counts.tsv"
done
echo "✓ Read-count table: $TRIM_DIR/trimmed_read_counts.tsv"

echo "All done."
echo "Raw FASTQs:     $RAW_DIR"
echo "Trimmed FASTQs: $TRIM_DIR"
echo "QC outputs:     $QC_DIR"
