#!/usr/bin/env bash
# SRA -> FASTQ pipeliine
# 1) prefetch .sra locally  2) vdb-validate  3) fasterq-dump (paired)  4) pigz/gzip  5) verify
# bifs619 group 5 
set -euo pipefail

# =======================
# CONFIG (override via env)
# =======================
THREADS="${THREADS:-8}"                      # CPUs for fasterq-dump & pigz
ROOT="${ROOT:-$HOME/sra_work}"               # base working directory
CACHE="${CACHE:-$ROOT/ncbi_sra_cache}"       # .sra cache
OUT="${OUT:-$ROOT/fastq}"                    # final FASTQs
TMP="${TMP:-$ROOT/tmp}"                      # temp dir (same filesystem as OUT recommended)
MAX_PREFETCH_GB="${MAX_PREFETCH_GB:-200}"    # per-run limit for prefetch
RETRIES="${RETRIES:-3}"                      # retry count
MINLEN="${MINLEN:-20}"                       # min read length to keep
PAIRED="--split-files"                       # paired-end (drop if single-end)

# =======================
# PREP
# =======================
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 SRR123456 SRR234567 ..."
  exit 1
fi

mkdir -p "$CACHE" "$OUT" "$TMP"

# point SRA Toolkit cache root (no-op if already set)
vdb-config -s "/repository/user/main/public/root=$CACHE" >/dev/null 2>&1 || true

echo "SRA Toolkit:"
command -v prefetch >/dev/null && prefetch --version || echo "prefetch not found"
command -v fasterq-dump >/dev/null && fasterq-dump --version || echo "fasterq-dump not found"
command -v vdb-validate >/dev/null && vdb-validate --version || echo "vdb-validate not found"

# generic retry helper
retry() {
  local max=$1; shift
  local n=1
  until "$@"; do
    if (( n >= max )); then
      echo "ERROR: command failed after $max attempts: $*"
      return 1
    fi
    echo "✗ Attempt $n failed. Retrying in $((n*5))s…" >&2
    sleep $((n*5))
    ((n++))
  done
}

# =======================
# MAIN
# =======================
for SRR in "$@"; do
  echo "=========================== $SRR ==========================="

  # --- 0) Normalize obvious typos like missing 'R' ---
  if [[ ! "$SRR" =~ ^SRR[0-9]+$ ]]; then
    echo "WARN: Accession '$SRR' doesn't look like SRR*; continuing but double-check it."
  fi

  # --- 1) Prefetch (local .sra) ---
  retry "$RETRIES" prefetch --progress -O "$CACHE" -X "${MAX_PREFETCH_GB}G" "$SRR"

  # resolve SRA path (sometimes sits at $CACHE/$SRR/$SRR.sra or $CACHE/$SRR.sra)
  SRA_PATH="$CACHE/$SRR/$SRR.sra"
  [[ -f "$SRA_PATH" ]] || SRA_PATH="$CACHE/$SRR.sra"
  if [[ ! -f "$SRA_PATH" ]]; then
    echo "ERROR: Could not find $SRR.sra under $CACHE after prefetch."
    exit 1
  fi

  # --- 2) Validate the container (.sra) ---
  echo "→ Validating $SRA_PATH"
  retry "$RETRIES" vdb-validate "$SRA_PATH"

  # --- 3) Idempotent guards (skip/compress/clean) ---
  if [[ -f "$OUT/${SRR}_1.fastq.gz" && -f "$OUT/${SRR}_2.fastq.gz" ]]; then
    echo "→ Outputs already exist (gz): skipping $SRR"
    continue
  fi

  if [[ -f "$OUT/${SRR}_1.fastq" && -f "$OUT/${SRR}_2.fastq" ]]; then
    echo "→ Found uncompressed FASTQs; compressing"
    if command -v pigz >/dev/null 2>&1; then
      pigz -p "$THREADS" "$OUT/${SRR}_1.fastq" "$OUT/${SRR}_2.fastq"
    else
      gzip -1 "$OUT/${SRR}_1.fastq" "$OUT/${SRR}_2.fastq"
    fi
    gzip -t "$OUT/${SRR}_1.fastq.gz" "$OUT/${SRR}_2.fastq.gz"
    echo "✓ $SRR compressed; skipping conversion"
    continue
  fi

  # remove any stray single-file artifact that blocks conversion
  rm -f "$OUT/${SRR}.fastq"

  # --- 4) Convert with fasterq-dump (paired, force overwrite) ---
  echo "→ Converting to FASTQ"
  retry "$RETRIES" fasterq-dump \
    -e "$THREADS" -p \
    -t "$TMP" \
    -O "$OUT" \
    $PAIRED --skip-technical --min-read-len "$MINLEN" \
    -f \
    "$SRA_PATH"

  # --- 5) Compress with pigz if available (fallback to gzip) ---
  echo "→ Compressing FASTQs"
  shopt -s nullglob
  FQS=( "$OUT/${SRR}_1.fastq" "$OUT/${SRR}_2.fastq" )
  if (( ${#FQS[@]} )); then
    if command -v pigz >/dev/null 2>&1; then
      pigz -p "$THREADS" "${FQS[@]}"
    else
      echo "  pigz not found; using gzip -1"
      gzip -1 "${FQS[@]}"
    fi
  fi

  # --- 6) Verify gzip integrity ---
  echo "→ Verifying gzip integrity"
  gzip -t "$OUT/${SRR}_1.fastq.gz" "$OUT/${SRR}_2.fastq.gz"

  # --- 7) Tiny summary (optional but handy) ---
  {
    r1_lines=$(zcat "$OUT/${SRR}_1.fastq.gz" | wc -l || echo 0)
    r2_lines=$(zcat "$OUT/${SRR}_2.fastq.gz" | wc -l || echo 0)
    r1=$(( r1_lines / 4 ))
    r2=$(( r2_lines / 4 ))
    printf "%s\t%s\t%s\n" "$SRR" "$r1" "$r2"
  } >> "$OUT/srr_read_counts.tsv"

  echo "✓ Done: $SRR"
done

echo "All requested SRRs completed. FASTQs in: $OUT"
echo "Per-SRR read counts appended to: $OUT/srr_read_counts.tsv"

