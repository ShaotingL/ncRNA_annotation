#!/usr/bin/env python3
"""
Convert Infernal cmscan --fmt 2 --tblout to GFF3.
cmscan fmt2 column layout (0-indexed):
  0:idx  1:query_name(CM)  2:query_acc  3:target(scaffold)  4:target_acc
  5:clan  6:mdl  7:mdl_from  8:mdl_to  9:seq_from  10:seq_to  11:strand
  12:trunc  13:pass  14:gc  15:bias  16:score  17:E-value  18:inc  19:olp
  20-25:window_info  26:mdl_len  27:seq_len  28+:description
Usage: python cmscan_tbl2gff3.py cmscan.tbl > output.gff3
"""
import sys

def guess_type(rfam_id, desc):
    d = (rfam_id + ' ' + desc).lower()
    if 'mir' in rfam_id.lower() or 'microrna' in d or 'mirna' in d:
        return 'miRNA'
    if 'snor' in d or 'snorna' in d:
        return 'snoRNA'
    if 'snrna' in d or 'spliceosomal' in d or rfam_id in ('U1','U2','U4','U5','U6','U11','U12','U4atac','U6atac'):
        return 'snRNA'
    if 'rrna' in d or 'ribosomal' in d or 'lsu_' in d or 'ssu_' in d or rfam_id.startswith('5S') or '5_8s' in d:
        return 'rRNA'
    if 'trna' in d:
        return 'tRNA'
    return 'ncRNA'

print("##gff-version 3")
counter = 0
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    parts = line.split()
    if len(parts) < 20:
        continue

    rfam_id  = parts[1]   # CM / query name
    rfam_acc = parts[2]   # CM accession (RF#####)
    target   = parts[3]   # genome scaffold (target sequence)
    clan     = parts[5]
    seq_from = int(parts[9])
    seq_to   = int(parts[10])
    strand   = parts[11]
    score    = parts[16]
    evalue   = parts[17]
    inc      = parts[18]  # ! = above GA threshold
    olp      = parts[19]  # x = excluded by clan competition

    # Skip clan-competition losers and below-threshold hits
    if olp == 'x':
        continue
    if inc not in ('!', '?'):
        continue

    start = min(seq_from, seq_to)
    end   = max(seq_from, seq_to)
    desc  = ' '.join(parts[28:]) if len(parts) > 28 else ''

    ncrna_type = guess_type(rfam_id, desc)
    counter += 1
    feat_id = f"{ncrna_type}_{counter:06d}"

    attrs = (f"ID={feat_id};Name={rfam_id};rfam_acc={rfam_acc};"
             f"clan={clan};score={score};evalue={evalue};"
             f"description={desc.replace(';','_')}")

    print('\t'.join([
        target, 'cmscan', ncrna_type,
        str(start), str(end), score, strand, '.', attrs
    ]))
