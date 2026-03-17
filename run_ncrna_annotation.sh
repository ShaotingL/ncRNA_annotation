#!/usr/bin/env bash
# ncRNA annotation pipeline — plant-specific Rfam subset
# Tools: tRNAscan-SE (tRNA) | cmscan+Rfam_plant+clanin (miRNA/snRNA/snoRNA) | barrnap (rRNA)
set -euo pipefail

CONDA_BASE=/home/st/miniconda3
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate bio-env

GENOME_DIR=/mnt/bay3_1.2T/2.GFAP/genome
RFAM_DIR=/mnt/bay3_1.2T/2.GFAP/rfam
OUT_BASE=/mnt/bay3_1.2T/2.GFAP/results
CPU=20
LOG="$OUT_BASE/ncrna_pipeline.log"

HAPA_GENOME="$GENOME_DIR/CR0040_HaplotypeA.fasta"
HAPB_GENOME="$GENOME_DIR/CR0040_HaplotypeB.fasta"
RFAM_CM="$RFAM_DIR/Rfam_plant.cm"        # 632 plant-relevant CMs only
RFAM_CLANIN="$RFAM_DIR/Rfam_plant.clanin"

exec &> >(tee -a "$LOG")
echo "========================================"
echo "ncRNA Annotation Pipeline (plant Rfam subset, 632 CMs)"
echo "Started: $(date)"
echo "========================================"

annotate_hap() {
    local HAP=$1
    local GENOME=$2
    local OUTDIR="$OUT_BASE/ncRNA_${HAP}"
    mkdir -p "$OUTDIR"

    echo ""
    echo "========================================"
    echo "  $HAP  |  $(date)"
    echo "========================================"

    # ── 1. tRNA — tRNAscan-SE ─────────────────────────────────────────────
    local TRNA_GFF="$OUTDIR/${HAP}_tRNA.gff3"
    if [[ ! -s "$TRNA_GFF" ]]; then
        echo "[$(date +%H:%M:%S)] [$HAP] tRNAscan-SE..."
        tRNAscan-SE -E --thread "$CPU" \
            -o "$OUTDIR/${HAP}_tRNA.txt" \
            --gff "$TRNA_GFF" \
            --log "$OUTDIR/${HAP}_tRNAscan.log" \
            "$GENOME"
        echo "[$(date +%H:%M:%S)] [$HAP] tRNA done: $(grep -vc '^#' "$TRNA_GFF") tRNAs"
    else
        echo "[$(date +%H:%M:%S)] [$HAP] tRNA exists ($(grep -vc '^#' "$TRNA_GFF")), skipping."
    fi

    # ── 2. rRNA — barrnap ─────────────────────────────────────────────────
    local RRNA_GFF="$OUTDIR/${HAP}_rRNA.gff3"
    if [[ ! -s "$RRNA_GFF" ]]; then
        echo "[$(date +%H:%M:%S)] [$HAP] barrnap (rRNA)..."
        # Prepend dummy ACGT sequence to ensure nhmmer can guess alphabet
        local GENOME_BARRNAP="/tmp/${HAP}_barrnap_input.fa"
        echo -e ">dummy_acgt_alphabet\nACGTACGTACGTACGTACGT" > "$GENOME_BARRNAP"
        cat "$GENOME" >> "$GENOME_BARRNAP"
        barrnap --kingdom euk --threads "$CPU" \
            --outseq "$OUTDIR/${HAP}_rRNA.fa" \
            "$GENOME_BARRNAP" 2>/dev/null | \
            grep -v "dummy_acgt" > "$RRNA_GFF"
        rm -f "$GENOME_BARRNAP"
        echo "[$(date +%H:%M:%S)] [$HAP] rRNA done: $(grep -vc '^#' "$RRNA_GFF") rRNAs"
    else
        echo "[$(date +%H:%M:%S)] [$HAP] rRNA exists, skipping."
    fi

    # ── 3. miRNA/snRNA/snoRNA — cmscan + plant Rfam subset ────────────────
    local CM_TBL="$OUTDIR/${HAP}_cmscan.tbl"
    local CM_GFF="$OUTDIR/${HAP}_rfam.gff3"
    if [[ ! -s "$CM_TBL" ]]; then
        echo "[$(date +%H:%M:%S)] [$HAP] cmscan (plant Rfam 632 CMs)..."
        cmscan \
            --cpu "$CPU" \
            --cut_ga \
            --rfam \
            --nohmmonly \
            --fmt 2 \
            --clanin "$RFAM_CLANIN" \
            --tblout "$CM_TBL" \
            -o "$OUTDIR/${HAP}_cmscan_full.out" \
            "$RFAM_CM" "$GENOME"
        echo "[$(date +%H:%M:%S)] [$HAP] cmscan done: $(grep -vc '^#' "$CM_TBL") hits"
    else
        echo "[$(date +%H:%M:%S)] [$HAP] cmscan exists, skipping."
    fi

    # ── 4. Convert to GFF3 ────────────────────────────────────────────────
    if [[ -s "$CM_TBL" && ! -s "$CM_GFF" ]]; then
        echo "[$(date +%H:%M:%S)] [$HAP] Converting to GFF3..."
        python3 /mnt/bay3_1.2T/2.GFAP/scripts/cmscan_tbl2gff3.py "$CM_TBL" > "$CM_GFF"
        echo "[$(date +%H:%M:%S)] [$HAP] GFF3: $(grep -vc '^#' "$CM_GFF") entries"
    fi

    echo ""
    echo "=== $HAP Summary ==="
    echo "  tRNA : $(grep -vc '^#' "$OUTDIR/${HAP}_tRNA.gff3"  2>/dev/null || echo '?')"
    echo "  rRNA : $(grep -vc '^#' "$OUTDIR/${HAP}_rRNA.gff3"  2>/dev/null || echo '?')"
    echo "  Rfam : $(grep -vc '^#' "$CM_GFF"                   2>/dev/null || echo 'pending')"
    echo "  Finished: $(date)"
}

annotate_hap "HapA" "$HAPA_GENOME"
annotate_hap "HapB" "$HAPB_GENOME"

echo ""
echo "========================================"
echo "All done: $(date)"
echo "========================================"
