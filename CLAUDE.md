# CLAUDE.md — ncRNA Annotation Project

This file documents the workflow for future reference and reproducibility.

## Project Overview

Genome-wide ncRNA annotation of *Vanilla planifolia* CR0040 (HaplotypeA & HaplotypeB)
using rigorous bioinformatics tools. Results are in `ncRNA_HapA/` and `ncRNA_HapB/`.

**Genome source**: CNCB TCOD database
- HaplotypeA: `https://download.cncb.ac.cn/tcod/reference/33/CR0040_HaplotypeA.tar.gz`
- HaplotypeB: `https://download.cncb.ac.cn/tcod/reference/34/CR0040_HaplotypeB.tar.gz`
- Genome FASTA stored at: `/mnt/bay3_1.2T/2.GFAP/genome/`

---

## Why This Method (vs GFAP)

| | GFAP | This workflow |
|--|--|--|
| miRNA model | Profile HMM (sequence only) | Covariance Model (sequence + **RNA secondary structure**) |
| Input | Transcripts (annotated only) | **Genome** (finds unannotated loci too) |
| Database | miRBase-derived | **Rfam** (manually curated, peer-reviewed) |
| rRNA | hmmsearch | **barrnap** (dedicated rRNA tool) |
| tRNA | hmmsearch | **tRNAscan-SE** (gold standard) |

Covariance models (Infernal/cmscan) are the community standard for ncRNA detection
because ncRNA function is defined by secondary structure, not sequence alone.

---

## Tools

| Tool | Version | Purpose | Conda env |
|------|---------|---------|-----------|
| tRNAscan-SE | 2.0.12 | tRNA annotation | bio-env |
| barrnap | 0.9 | rRNA annotation (uses nhmmer internally) | bio-env |
| Infernal cmscan | 1.1.5 | miRNA/snRNA/snoRNA via Rfam CMs | bio-env |
| hmmer (nhmmer) | 3.4 | Dependency of barrnap | bio-env |

Install all tools:
```bash
conda install -n bio-env -c bioconda infernal trnascan-se barrnap hmmer -y
chmod +x ~/miniconda3/envs/bio-env/bin/nhmmer
chmod +x ~/miniconda3/envs/bio-env/bin/hmmpress
chmod +x ~/miniconda3/envs/bio-env/bin/makehmmerdb
# Fix barrnap's nhmmer lookup path
mkdir -p ~/miniconda3/envs/bio-env/lib/barrnap/binaries/linux
ln -sf ~/miniconda3/envs/bio-env/bin/nhmmer \
       ~/miniconda3/envs/bio-env/lib/barrnap/binaries/linux/nhmmer
```

---

## Rfam Database Setup

```bash
cd /mnt/bay3_1.2T/3.ncRNA_annotation/rfam/

# Download full Rfam CM database (~315 MB)
wget -c "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" -O Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm

# Download clan file (for clan competition / removing redundant hits)
wget "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin" -O Rfam.clanin
```

### Why a Plant-Specific Subset?

Full Rfam has 4,227 CMs, most are bacterial/metazoan-specific. Scanning a 1.4–2.0 Gb
plant genome against all 4,227 CMs takes weeks. We extracted 632 plant-relevant CMs:

| Category | Count |
|----------|-------|
| Plant miRNA (MIR*) | 313 |
| snoRNA | 271 |
| snRNA + plant-specific | 48 |
| **Total** | **632** |

```bash
# Regenerate plant subset (if needed)
source ~/miniconda3/etc/profile.d/conda.sh && conda activate bio-env

cmstat Rfam.cm | grep -v "^#" | awk '{print $2}' | python3 -c "
import sys
keep = []
for n in sys.stdin.read().split():
    nl = n.lower()
    if n.startswith('MIR'): keep.append(n)
    elif nl.startswith('snor') or nl.startswith('snod'): keep.append(n)
    elif n in ('U1','U2','U4','U5','U6','U11','U12','U4atac','U6atac','U3','U54','RNase_MRP'): keep.append(n)
    elif any(x in nl for x in ['plant_srp','plant_u','plastid','rnase_mrp','acea_u3']): keep.append(n)
for n in keep: print(n)
" > plant_cm_list.txt

cmfetch -f Rfam.cm plant_cm_list.txt > Rfam_plant.cm
cmpress Rfam_plant.cm

# Filtered clanin (only clans whose members are in our subset)
python3 -c "
plant = set(open('plant_cm_list.txt').read().split())
out = []
for line in open('Rfam.clanin'):
    parts = line.strip().split()
    members = [m for m in parts[1:] if m in plant]
    if members: out.append(parts[0] + '\t' + '\t'.join(members))
open('Rfam_plant.clanin','w').write('\n'.join(out)+'\n')
"
```

---

## Workflow

Run the full pipeline (HapA then HapB automatically):

```bash
screen -S ncrna bash -c 'bash /mnt/bay3_1.2T/3.ncRNA_annotation/run_ncrna_annotation.sh 2>&1'
# Detach: Ctrl+A D   |   Reattach: screen -r ncrna
```

The script skips completed steps automatically (safe to resume after interruption).

### Step 1 — tRNA (tRNAscan-SE)

```bash
tRNAscan-SE -E --thread 20 \
    -o ${HAP}_tRNA.txt \
    --gff ${HAP}_tRNA.gff3 \
    --log ${HAP}_tRNAscan.log \
    genome.fasta
```

- `-E`: eukaryotic mode
- Uses Infernal CM internally for high accuracy

### Step 2 — rRNA (barrnap)

```bash
# NOTE: barrnap fails on sequences with no G in the first few kb (e.g. telomeric scaffolds).
# Workaround: prepend a short ACGT dummy sequence so nhmmer can detect the alphabet.
echo -e ">dummy_acgt\nACGTACGTACGTACGTACGT" > input_with_dummy.fa
cat genome.fasta >> input_with_dummy.fa

barrnap --kingdom euk --threads 20 \
    --outseq ${HAP}_rRNA.fa \
    input_with_dummy.fa 2>/dev/null | \
    grep -v "dummy_acgt" > ${HAP}_rRNA.gff3
```

- `--kingdom euk`: eukaryotic rRNA models (18S, 5.8S, 28S, 5S)

### Step 3 — miRNA / snRNA / snoRNA (cmscan + Rfam)

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
- `--cut_ga`: use Rfam gathering thresholds (family-specific, peer-reviewed cutoffs)
- `--rfam`: activate genome scanning mode with HMM filter layer
- `--nohmmonly`: require CM evaluation (no HMM-only hits)
- `--fmt 2`: extended tblout format with clan competition columns
- `--clanin`: clan file to remove redundant hits from the same RNA family clan

### Step 4 — Convert to GFF3

```bash
python3 cmscan_tbl2gff3.py ${HAP}_cmscan.tbl > ${HAP}_rfam.gff3
```

The script reads fmt2 tblout and:
- Skips clan competition losers (`olp == 'x'`)
- Skips below-threshold hits (`inc` not in `!` or `?`)
- Infers ncRNA type from family name
- Outputs standard GFF3

### Extract miRNA only (from full rfam results)

```bash
grep $'\tmiRNA\t' ncRNA_HapA/HapA_rfam.gff3 > ncRNA_HapA/HapA_miRNA.gff3
grep $'\tmiRNA\t' ncRNA_HapB/HapB_rfam.gff3 > ncRNA_HapB/HapB_miRNA.gff3
```

---

## miRNA-Only Workflow (Faster)

If you only need miRNA annotation, use the dedicated script with **313 CMs** (miRNA only).
This is ~2× faster than the full plant workflow (632 CMs) and skips tRNA/rRNA steps.

### Rfam subset
- `rfam/Rfam_mirna.cm` — 313 plant miRNA CMs (MIR* families only), 25 MB
- `rfam/Rfam_mirna.clanin` — filtered clan file

### Run

```bash
screen -S mirna bash -c '
bash /mnt/bay3_1.2T/3.ncRNA_annotation/run_mirna_only.sh \
    /mnt/bay3_1.2T/2.GFAP/genome/CR0040_HaplotypeA.fasta HapA 20
bash /mnt/bay3_1.2T/3.ncRNA_annotation/run_mirna_only.sh \
    /mnt/bay3_1.2T/2.GFAP/genome/CR0040_HaplotypeB.fasta HapB 20
'
```

Results saved to `miRNA_results/`:
- `{PREFIX}_miRNA.gff3` — genomic coordinates of miRNA loci
- `{PREFIX}_cmscan.tbl` — raw cmscan output (fmt2)

### Apply to any plant genome

```bash
bash run_mirna_only.sh <genome.fasta> <prefix> [cpu]
# Example for a new species:
bash run_mirna_only.sh /path/to/species.fasta MySpecies 20
```

### Speed estimate (20 CPU)

| Genome size | Estimated time |
|-------------|---------------|
| ~500 Mb | ~1–2 hours |
| ~1.5 Gb | ~3–4 hours |
| ~2.0 Gb | ~4–6 hours |

> Note: Large repetitive scaffolds (unanchored contigs) take disproportionately longer
> due to high HMM filter hit rates from repetitive sequences.

---

## Results Summary

| Type | Tool | HapA | HapB |
|------|------|------|------|
| tRNA | tRNAscan-SE | 1,436 | 4,002 |
| rRNA | barrnap | 7,398 | 9,253 |
| miRNA | cmscan (Rfam_plant) | 87 | 100 |
| snRNA/snoRNA | cmscan (Rfam_plant) | 348 | 806 |

HapB has more features than HapA due to larger genome size (1.97 Gb vs 1.42 Gb)
and a large repetitive scaffold (VANPL_B_00000).

---

## Known Issues & Workarounds

### barrnap fails on HapB (telomeric scaffold)
VANPL_B_00001 starts with pure CTAATCCC telomeric repeats (no G in first ~1 kb).
nhmmer cannot guess the DNA alphabet from this sequence.
**Workaround**: prepend a dummy `ACGT` sequence before running barrnap (see Step 2 above).

### cmscan very slow on large repetitive scaffolds
VANPL_A_00000 (876 Mb) and VANPL_B_00000 are large unanchored scaffolds with
high repeat content, causing excessive HMM filter hits and very slow scanning.
With the 632-CM plant subset, each ~90 Mb chunk takes ~9 hours.
Full Rfam (4,227 CMs) would be ~6× slower and is not recommended for this genome.

### nhmmer / hmmpress permission errors
Conda sometimes installs HMMER binaries without execute permission.
```bash
chmod +x ~/miniconda3/envs/bio-env/bin/nhmmer
chmod +x ~/miniconda3/envs/bio-env/bin/hmmpress
```

---

## Reference

- Infernal/cmscan: Nawrocki & Eddy, *Bioinformatics* 2013. doi:10.1093/bioinformatics/btt509
- Rfam: Kalvari et al., *Nucleic Acids Research* 2021. doi:10.1093/nar/gkaa1047
- tRNAscan-SE: Chan et al., *Nucleic Acids Research* 2021. doi:10.1093/nar/gkab688
- barrnap: Seemann T. https://github.com/tseemann/barrnap
