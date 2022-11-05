from pathlib import Path

import numpy as np
from Bio import SeqIO

# Path to aligned fasta file
FILE_NAME = Path('../GitHub/sars2oolkit/sarbeco_family_aligned.fasta')
LIMIT = 23
REF_GENOME = "NC_045512.2"
np.seterr(divide='ignore', invalid='ignore')


def ungapped_pos(seq, pos):
    # Counter for gaps (-) and non gaps (a c g t)
    non_gap = 0
    gaps = 0
    # Loop over each nucleotide in the sequence
    for nt in seq:
        # If nucleotide is not a gap add 1 to non_gap
        if nt != '-':
            non_gap += 1
        # Else add 1 to gaps
        else:
            gaps += 1
        # Return position without gaps
        if pos == (gaps + non_gap):
            return pos - gaps


# Create a list of SeqIO Records
fasta = list(SeqIO.parse(FILE_NAME, format='fasta'))


# Create a dictionary to access each genome using 'name'
seqs = {}
for entry in fasta:
    seqs[entry.id] = entry


def conserved():
    position_legend = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    count_a = 0
    count_c = 0
    count_g = 0
    count_t = 0
    nucleotide = 1
    conserved_list = []
    while nucleotide <= len(seqs[REF_GENOME].seq):
        for record in fasta:
            nt = record.seq[nucleotide - 1]
            if nt != '-':
                if nt == 'a':
                    count_a += 1
                elif nt == 'c':
                    count_c += 1
                elif nt == 'g':
                    count_g += 1
                elif nt == 't':
                    count_t += 1
        data = np.array([count_a, count_c, count_g, count_t])
        percentages = data / data.sum(axis=0)
        out = list(zip(data, percentages))

        for i, item in enumerate(out):
            reference = seqs[REF_GENOME].seq[nucleotide - 1]
            if item[0] >= LIMIT and reference != '-':
                check_ref_nucleotide = position_legend[reference]
                if check_ref_nucleotide != i:
                    conserved_list.append(
                        [int(ungapped_pos(seqs[REF_GENOME].seq, nucleotide)), i])
        count_a = count_c = count_g = count_t = 0
        nucleotide += 1

    return conserved_list
