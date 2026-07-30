# amr_pure.py - Protein Feature Engineering Pipeline
# Pure Python implementation (No pip / No external packages)

import os
import csv
from itertools import product
from collections import defaultdict

print("[1/5] Initializing FASTA parser for Protein Sequences...")

def parse_fasta(filepath):
    sequences = []
    current_header = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_header:
            sequences.append((current_header, ''.join(current_seq)))
            
    return sequences

# Kyte-Doolittle Hydropathy Scale for GRAVY score calculation
HYDROPATHY_SCALE = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")
DIPEPTIDES = [''.join(p) for p in product(STANDARD_AA, repeat=2)]

fasta_files = {
    'KPC': 'bla_KPC.fasta',
    'NDM': 'bla_NDM.fasta',
    'TEM': 'bla_TEM.fasta',
    'OXA': 'bla_OXA.fasta'
}

data = []

for family, path in fasta_files.items():
    if not os.path.exists(path):
        continue
    records = parse_fasta(path)
    for seq_id, seq in records:
        # Filter valid protein sequences
        clean_seq = ''.join([aa for aa in seq if aa in STANDARD_AA])
        if len(clean_seq) >= 30:
            data.append({
                'ID': seq_id,
                'Family': family,
                'Length': len(clean_seq),
                'Sequence': clean_seq
            })

print(f"Loaded {len(data)} valid protein sequences.")

print("[2/5] Calculating Physicochemical Properties...")

for row in data:
    seq = row['Sequence']
    length = row['Length']
    
    # 1. GRAVY Score (Hydrophobicity)
    gravy = sum(HYDROPATHY_SCALE.get(aa, 0) for aa in seq) / length
    
    # 2. Net Charge at pH 7.0
    pos_charge = seq.count('K') + seq.count('R') + (0.1 * seq.count('H'))
    neg_charge = seq.count('D') + seq.count('E')
    net_charge = pos_charge - neg_charge
    
    # 3. Basic vs Acidic Ratio
    basic_ratio = (seq.count('K') + seq.count('R') + seq.count('H')) / length
    acidic_ratio = (seq.count('D') + seq.count('E')) / length
    
    # 4. Aromaticity (F, W, Y)
    aromaticity = (seq.count('F') + seq.count('W') + seq.count('Y')) / length
    
    row['GRAVY'] = round(gravy, 4)
    row['NetCharge_pH7'] = round(net_charge, 2)
    row['Basic_Ratio'] = round(basic_ratio, 4)
    row['Acidic_Ratio'] = round(acidic_ratio, 4)
    row['Aromaticity'] = round(aromaticity, 4)

print("[3/5] Computing 400 Dipeptide Frequencies (20x20 Matrix)...")

def compute_dipeptides(seq):
    total = len(seq) - 1
    if total <= 0:
        return {dp: 0.0 for dp in DIPEPTIDES}
    counts = defaultdict(int)
    for i in range(total):
        dp = seq[i:i+2]
        if dp in counts or dp in DIPEPTIDES:
            counts[dp] += 1
    return {dp: round(counts[dp] / total, 6) for dp in DIPEPTIDES}

for row in data:
    dp_freqs = compute_dipeptides(row['Sequence'])
    row.update(dp_freqs)

print("[4/5] Writing Dataset to CSV...")

metadata_cols = ['ID', 'Family', 'Length', 'GRAVY', 'NetCharge_pH7', 'Basic_Ratio', 'Acidic_Ratio', 'Aromaticity']
all_headers = metadata_cols + DIPEPTIDES

with open('amr_protein_features.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=all_headers, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(data)

print("[5/5] Generating Summary Report...\n")

family_summary = defaultdict(lambda: {
    'count': 0, 'len': 0, 'gravy': 0, 'charge': 0, 'aromatic': 0
})

for d in data:
    fam = d['Family']
    family_summary[fam]['count'] += 1
    family_summary[fam]['len'] += d['Length']
    family_summary[fam]['gravy'] += d['GRAVY']
    family_summary[fam]['charge'] += d['NetCharge_pH7']
    family_summary[fam]['aromatic'] += d['Aromaticity']

print("=" * 75)
print(f"{'FAMILY':<8} | {'COUNT':<6} | {'AVG LEN':<8} | {'AVG GRAVY':<10} | {'NET CHARGE':<11} | {'AROMATIC%':<9}")
print("=" * 75)
for fam, s in family_summary.items():
    c = s['count']
    print(f"{fam:<8} | {c:<6} | {round(s['len']/c,1):<8} | {round(s['gravy']/c,3):<10} | {round(s['charge']/c,2):<11} | {round((s['aromatic']/c)*100,2)}%")
print("=" * 75)

print("\nProcess Complete!")
print("Output generated: 'amr_protein_features.csv' (Contains 408 total feature columns).")