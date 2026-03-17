#!/usr/bin/env bash
# miRNA-only annotation pipeline using Infernal cmscan + Rfam plant miRNA CMs (313 CMs)
# Usage: bash run_mirna_only.sh <genome.fasta> <output_prefix> [cpu]
# Example:
#   bash run_mirna_only.sh /mnt/bay3_1.2T/2.GFAP/genome/CR0040_HaplotypeA.fasta HapA 20
#   bash run_mirna_only.sh /mnt/bay3_1.2T/2.GFAP/genome/CR0040_HaplotypeB.fasta HapB 20
set -euo pipefail

CONDA_BASE=/home/st/miniconda3
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate bio-env

SCRIPT_DIR=/mnt/bay3_1.2T/3.ncRNA_annotation
RFAM_CM="$SCRIPT_DIR/rfam/Rfam_mirna.cm"
RFAM_CLANIN="$SCRIPT_DIR/rfam/Rfam_mirna.clanin"
OUTDIR="$SCRIPT_DIR/miRNA_results"

GENOME=${1:?Usage: bash run_mirna_only.sh <genome.fasta> <prefix> [cpu]}
PREFIX=${2:?Usage: bash run_mirna_only.sh <genome.fasta> <prefix> [cpu]}
CPU=${3:-20}

mkdir -p "$OUTDIR"
LOG="$OUTDIR/${PREFIX}_mirna.log"

exec &> >(tee -a "$LOG")
echo "========================================"
echo "miRNA Annotation — Rfam plant subset (313 CMs)"
echo "  Genome : $GENOME"
echo "  Prefix : $PREFIX"
echo "  CPU    : $CPU"
echo "  Started: $(date)"
echo "========================================"

TBL="$OUTDIR/${PREFIX}_cmscan.tbl"
GFF="$OUTDIR/${PREFIX}_miRNA.gff3"

# ── cmscan ────────────────────────────────────────────────────────────────────
if [[ ! -s "$TBL" ]]; then
    echo "[$(date +%H:%M:%S)] Running cmscan (313 miRNA CMs)..."
    cmscan \
        --cpu "$CPU" \
        --cut_ga \
        --rfam \
        --nohmmonly \
        --fmt 2 \
        --clanin "$RFAM_CLANIN" \
        --tblout "$TBL" \
        -o "$OUTDIR/${PREFIX}_cmscan_full.out" \
        "$RFAM_CM" "$GENOME"
    echo "[$(date +%H:%M:%S)] cmscan done: $(grep -vc '^#' "$TBL") hits"
else
    echo "[$(date +%H:%M:%S)] tblout exists ($(grep -vc '^#' "$TBL") hits), skipping."
fi

# ── Convert to GFF3 ───────────────────────────────────────────────────────────
echo "[$(date +%H:%M:%S)] Converting to GFF3..."
python3 "$SCRIPT_DIR/cmscan_tbl2gff3.py" "$TBL" > "$GFF"
TOTAL=$(grep -vc '^#' "$GFF")
echo "[$(date +%H:%M:%S)] GFF3 written: $TOTAL miRNA loci"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== $PREFIX miRNA Summary ==="
echo "  Total loci : $TOTAL"
echo "  Families   :"
grep -v '^#' "$GFF" | grep -o 'Name=[^;]*' | sort | uniq -c | sort -rn | \
    awk '{printf "    %-20s %s\n", $2, $1}'
echo "  Output GFF3: $GFF"
echo "  Finished   : $(date)"
echo "========================================"
