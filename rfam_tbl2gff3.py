#!/usr/bin/env python3
"""
Convert Infernal cmsearch --tblout to GFF3.
Usage: python rfam_tbl2gff3.py cmsearch.tbl > output.gff3
"""
import sys, re

# Rfam family → ncRNA type mapping (common families)
TYPE_MAP = {
    'mir': 'miRNA', 'miRNA': 'miRNA',
    'sno': 'snoRNA', 'snRNA': 'snRNA', 'sRNA': 'snRNA',
    'rRNA': 'rRNA', 'LSU': 'rRNA', 'SSU': 'rRNA', '5S': 'rRNA', '5_8S': 'rRNA',
    'tRNA': 'tRNA',
}

def guess_type(rfam_id, rfam_desc):
    desc_lower = rfam_desc.lower()
    if 'mirna' in desc_lower or 'microrna' in desc_lower or 'mir-' in desc_lower:
        return 'miRNA'
    if 'snorna' in desc_lower or 'sno' in desc_lower:
        return 'snoRNA'
    if 'snrna' in desc_lower or 'spliceosomal' in desc_lower:
        return 'snRNA'
    if 'rrna' in desc_lower or 'ribosomal' in desc_lower or 'lsu' in desc_lower or 'ssu' in desc_lower:
        return 'rRNA'
    if 'trna' in desc_lower:
        return 'tRNA'
    return 'ncRNA'

print("##gff-version 3")
counter = 0
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    # cmsearch --tblout columns (space-delimited, last field is description)
    # target_name, accession, query_name, accession, mdl, mdl_from, mdl_to,
    # seq_from, seq_to, strand, trunc, pass, gc, bias, score, E-value, inc, description
    parts = line.split()
    if len(parts) < 18:
        continue
    target   = parts[0]   # chromosome/scaffold
    rfam_acc = parts[2]   # RF#####
    rfam_id  = parts[3]   # Rfam family ID
    seq_from = int(parts[7])
    seq_to   = int(parts[8])
    strand   = parts[9]   # + or -
    score    = parts[14]
    evalue   = parts[15]
    desc     = ' '.join(parts[17:])

    if strand == '-':
        start, end = min(seq_from, seq_to), max(seq_from, seq_to)
    else:
        start, end = seq_from, seq_to

    ncrna_type = guess_type(rfam_acc, desc)
    counter += 1
    feat_id = f"{ncrna_type}_{counter:06d}"

    attrs = (f"ID={feat_id};Name={rfam_id};rfam_acc={rfam_acc};"
             f"score={score};evalue={evalue};description={desc.replace(';','_')}")

    print('\t'.join([
        target, 'cmsearch', ncrna_type,
        str(start), str(end), score, strand, '.', attrs
    ]))
