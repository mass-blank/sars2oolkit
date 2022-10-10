from functools import cache
from Bio import SeqIO
import numpy as np

from pathlib import Path

# Path to aligned fasta file
FILE_NAME = Path("/Users/nazar/SARS-CoV-2/sarbeco_family_aligned.fasta")
LIMIT = 16
REF_GENOME = "NC_045512.2"
np.seterr(divide='ignore', invalid='ignore')

# Check ungapped nucleotide position
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
# Create a list of SeqIO  
fasta = list(SeqIO.parse(FILE_NAME, format='fasta'))

# Create a dictionary to access each genome using 'name'
seqs = {}
for entry in fasta:
    seqs[entry.id] = entry

# numerate nucleotides
def conserved():
    position_legend = {'a': 0, 'c': 1, 'g': 2, 't': 3}
    countA = 0
    countC = 0 
    countG = 0
    countT = 0
    nucleotide = 1
    conserved_list = []
    with open('conserved.txt', 'w') as file:
        while nucleotide <= len(seqs[REF_GENOME].seq):
            for entry in fasta:
                nt = entry.seq[nucleotide-1]
                if nt != '-':
                    if nt == 'a':
                        countA += 1
                    elif nt == 'c':
                        countC += 1
                    elif nt == 'g':
                        countG += 1
                    elif nt == 't':
                        countT += 1
            data = np.array([countA, countC, countG, countT])
            percentages = data/data.sum(axis=0)
            out = list(zip(data, percentages))
                    
            for i, item in enumerate(out):
                reference = seqs[REF_GENOME].seq[nucleotide-1]
                if item[0] >= LIMIT and reference != '-':
                    check_ref_nucleotide = position_legend[reference]
                    if check_ref_nucleotide != i:
                        conserved_list.append([int(ungapped_pos(seqs[REF_GENOME].seq, nucleotide)), i])
                        # s = str(f"POS:\t{ungapped_pos(seqs[REF_GENOME].seq, nucleotide)}\tconserved\tREF: {str.capitalize(reference)}\tU_GAP:{nucleotide}\n")
                        # print(s)

            countA = countC = countG = countT = 0
            nucleotide += 1

    return conserved_list