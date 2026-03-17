# ncRNA Annotation — *Vanilla planifolia* CR0040

Genome-wide non-coding RNA (ncRNA) annotation of *Vanilla planifolia* CR0040 (HaplotypeA & HaplotypeB) using community-standard bioinformatics tools: **tRNAscan-SE**, **barrnap**, and **Infernal cmscan** with Rfam covariance models.

---

## Why This Approach?

ncRNA function is defined by **RNA secondary structure**, not sequence alone. This workflow uses **covariance models (CMs)** that capture both sequence conservation and structural constraints — the community standard for ncRNA detection.

| Feature | GFAP | This Workflow |
|---------|------|---------------|
| miRNA model | Profile HMM (sequence only) | **Covariance Model** (sequence + RNA secondary structure) |
| Input | Annotated transcripts only | **Genome** (finds unannotated loci too) |
| Database | miRBase-derived | **Rfam** (manually curated, peer-reviewed) |
| rRNA tool | hmmsearch | **barrnap** (dedicated rRNA tool) |
| tRNA tool | hmmsearch | **tRNAscan-SE** (gold standard) |

---

## Tools

| Tool | Version | Purpose | Conda env |
|------|---------|---------|-----------|
| tRNAscan-SE | 2.0.12 | tRNA annotation | `bio-env` |
| barrnap | 0.9 | rRNA annotation (18S, 5.8S, 28S, 5S) | `bio-env` |
| Infernal cmscan | 1.1.5 | miRNA / snRNA / snoRNA via Rfam CMs | `bio-env` |
| hmmer (nhmmer) | 3.4 | Dependency of barrnap | `bio-env` |

### Installation

```bash
conda install -n bio-env -c bioconda infernal trnascan-se barrnap hmmer -y

# Fix barrnap's nhmmer lookup path
mkdir -p ~/miniconda3/envs/bio-env/lib/barrnap/binaries/linux
ln -sf ~/miniconda3/envs/bio-env/bin/nhmmer \
       ~/miniconda3/envs/bio-env/lib/barrnap/binaries/linux/nhmmer

# Fix execute permissions (conda sometimes omits them)
chmod +x ~/miniconda3/envs/bio-env/bin/nhmmer
chmod +x ~/miniconda3/envs/bio-env/bin/hmmpress
```

---

## Repository Structure

```
ncRNA_annotation/
├── run_ncrna_annotation.sh      # Full pipeline (tRNA + rRNA + miRNA/snRNA/snoRNA)
├── run_mirna_only.sh            # miRNA-only pipeline (~2x faster)
├── cmscan_tbl2gff3.py           # Convert Infernal cmscan tblout -> GFF3
├── rfam/
│   ├── Rfam_plant.cm            # Plant-relevant Rfam subset (632 CMs)
│   ├── Rfam_plant.clanin        # Filtered clan file for plant subset
│   ├── Rfam_mirna.cm            # Plant miRNA only (313 CMs, MIR* families)
│   └── Rfam_mirna.clanin        # Filtered clan file for miRNA subset
│   # (full Rfam.cm not tracked — too large; see Setup below)
├── ncRNA_HapA/
│   ├── HapA_tRNA.gff3           # tRNA loci
│   ├── HapA_rRNA.gff3           # rRNA loci
│   ├── HapA_rfam.gff3           # miRNA + snRNA/snoRNA loci
│   └── HapA_cmscan.tbl          # Raw cmscan fmt2 tblout
├── ncRNA_HapB/
│   ├── HapB_tRNA.gff3
│   ├── HapB_rRNA.gff3
│   ├── HapB_rfam.gff3
│   └── HapB_cmscan.tbl
└── CLAUDE.md                    # Detailed protocol and troubleshooting notes
```

---

## Rfam Database Setup

The Rfam CM databases are not tracked in git (too large). Download and prepare them once:

```bash
cd rfam/

# Full Rfam (~315 MB) — only needed to regenerate plant subsets
wget -c "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" -O Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm

wget "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin" -O Rfam.clanin
```

The plant-specific subsets (`Rfam_plant.cm`, `Rfam_mirna.cm`) are already committed and pressed — no further setup needed for normal runs.

### Why a Plant Subset?

Full Rfam has 4,227 CMs, most bacterial/metazoan. Scanning a 1.4–2.0 Gb plant genome against all CMs would take weeks. We extracted plant-relevant families:

| Category | CMs |
|----------|-----|
| Plant miRNA (MIR\*) | 313 |
| snoRNA | 271 |
| snRNA + plant-specific | 48 |
| **Total (plant subset)** | **632** |

---

## Quick Start

### Full Annotation (tRNA + rRNA + miRNA/snRNA/snoRNA)

```bash
# Run in screen (takes several hours for Gb-scale genomes)
screen -S ncrna bash -c 'bash run_ncrna_annotation.sh 2>&1'
# Detach: Ctrl+A D   |   Reattach: screen -r ncrna
```

Edit the genome paths at the top of `run_ncrna_annotation.sh` to point to your FASTA files. The script runs HapA then HapB automatically and **skips completed steps** (safe to resume after interruption).

### miRNA-Only Annotation (Faster)

Uses 313 miRNA CMs only (~2x faster than the full plant workflow):

```bash
bash run_mirna_only.sh <genome.fasta> <output_prefix> [cpu]

# Examples:
bash run_mirna_only.sh /path/to/HaplotypeA.fasta HapA 20
bash run_mirna_only.sh /path/to/HaplotypeB.fasta HapB 20
```

Results are saved to `miRNA_results/`:
- `{PREFIX}_miRNA.gff3` — genomic coordinates
- `{PREFIX}_cmscan.tbl` — raw cmscan output
- `{PREFIX}_mirna.log` — run log

### Apply to Any Plant Genome

```bash
bash run_mirna_only.sh /path/to/species.fasta MySpecies 20
```

#### Speed Estimate (20 CPUs)

| Genome Size | Estimated Time |
|-------------|---------------|
| ~500 Mb | ~1–2 hours |
| ~1.5 Gb | ~3–4 hours |
| ~2.0 Gb | ~4–6 hours |

> Note: Large repetitive/unanchored scaffolds take disproportionately longer due to high HMM filter hit rates.

---

## Workflow Details

### Step 1 — tRNA (tRNAscan-SE)

```bash
tRNAscan-SE -E --thread 20 \
    -o ${HAP}_tRNA.txt \
    --gff ${HAP}_tRNA.gff3 \
    --log ${HAP}_tRNAscan.log \
    genome.fasta
```

`-E`: eukaryotic mode; uses Infernal CM internally for maximum sensitivity.

### Step 2 — rRNA (barrnap)

```bash
# Prepend dummy sequence to avoid nhmmer alphabet detection failure
# (needed for genomes starting with telomeric repeats lacking G)
echo -e ">dummy_acgt\nACGTACGTACGTACGTACGT" > input_with_dummy.fa
cat genome.fasta >> input_with_dummy.fa

barrnap --kingdom euk --threads 20 \
    --outseq ${HAP}_rRNA.fa \
    input_with_dummy.fa 2>/dev/null | \
    grep -v "dummy_acgt" > ${HAP}_rRNA.gff3
```

### Step 3 — miRNA/snRNA/snoRNA (Infernal cmscan + Rfam)

```bash
cmscan \
    --cpu 20 \
    --cut_ga \
    --rfam \
    --nohmmonly \
    --fmt 2 \
    --clanin rfam/Rfam_plant.clanin \
    --tblout ${HAP}_cmscan.tbl \
    -o ${HAP}_cmscan_full.out \
    rfam/Rfam_plant.cm \
    genome.fasta
```

Key flags:
- `--cut_ga` — Rfam gathering thresholds (peer-reviewed, family-specific cutoffs)
- `--rfam` — genome scanning mode with HMM pre-filter
- `--nohmmonly` — require CM-level evaluation (no HMM-only hits)
- `--fmt 2` — extended tblout with clan competition columns
- `--clanin` — remove redundant hits within the same RNA family clan

### Step 4 — Convert to GFF3

```bash
python3 cmscan_tbl2gff3.py ${HAP}_cmscan.tbl > ${HAP}_rfam.gff3
```

Filters out clan competition losers (`olp == 'x'`) and below-threshold hits.

---

## Results — *Vanilla planifolia* CR0040

| Type | Tool | HapA | HapB |
|------|------|------|------|
| tRNA | tRNAscan-SE | 1,436 | 4,002 |
| rRNA | barrnap | 7,398 | 9,253 |
| miRNA | cmscan (Rfam_plant) | 87 | 100 |
| snRNA/snoRNA | cmscan (Rfam_plant) | 348 | 806 |

HapB has more features than HapA due to larger genome size (1.97 Gb vs 1.42 Gb) and a large repetitive scaffold (VANPL_B_00000).

**Genome source**: CNCB TCOD database
- HaplotypeA: `https://download.cncb.ac.cn/tcod/reference/33/CR0040_HaplotypeA.tar.gz`
- HaplotypeB: `https://download.cncb.ac.cn/tcod/reference/34/CR0040_HaplotypeB.tar.gz`

---

## References

- **Infernal/cmscan**: Nawrocki & Eddy, *Bioinformatics* 2013. doi:[10.1093/bioinformatics/btt509](https://doi.org/10.1093/bioinformatics/btt509)
- **Rfam**: Kalvari et al., *Nucleic Acids Research* 2021. doi:[10.1093/nar/gkaa1047](https://doi.org/10.1093/nar/gkaa1047)
- **tRNAscan-SE**: Chan et al., *Nucleic Acids Research* 2021. doi:[10.1093/nar/gkab688](https://doi.org/10.1093/nar/gkab688)
- **barrnap**: Seemann T. [github.com/tseemann/barrnap](https://github.com/tseemann/barrnap)
