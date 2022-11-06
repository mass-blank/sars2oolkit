from pathlib import Path
import os
import numpy as np
from Bio import SeqIO

project_root = Path(__file__).parent.resolve()
# Path to aligned fasta file
FILE_NAME = Path(project_root,'sarbeco_family_aligned.fasta')
# Threshold for minimum number that each nucleotide position
# per genome must contain to be accepted as 'ancestral'
LIMIT = 23
REF_GENOME = "NC_045512.2"
np.seterr(divide='ignore', invalid='ignore')


def ungapped_pos(seq, pos):
    """This function returns the un-gapped position"""
    # Counter for gaps (-) and non-gaps (a c g t)
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


# Create a list of SeqIO Record objects from aligned .fasta file
fasta = list(SeqIO.parse(FILE_NAME, format='fasta'))


# Create a dictionary to access each genome using 'name' i.e. "NC_045512.2"
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
    # loop through each nucleotide position
    while nucleotide <= len(seqs[REF_GENOME].seq):
        # loop through each nucleotide per genome
        for record in fasta:
            # return nucleotide value at this position
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
        # create a numpy array with sum total per nucleotide
        data = np.array([count_a, count_c, count_g, count_t])
        # convert to float (proportion rather than absolute)
        percentages = data / data.sum(axis=0)
        # combine both into 1 data structure
        out = list(zip(data, percentages))

        for i, item in enumerate(out):
            # reference nucleotide value
            reference = seqs[REF_GENOME].seq[nucleotide - 1]
            # must be higher or equal to threshold and reference must not be a gapped position
            if item[0] >= LIMIT and reference != '-':
                check_ref_nucleotide = position_legend[reference]
                # if reference nucleotide does not match conserved position append that position and conserved allele
                # to a list
                if check_ref_nucleotide != i:
                    conserved_list.append(
                        [int(ungapped_pos(seqs[REF_GENOME].seq, nucleotide)), i])
        # clear counter
        count_a = count_c = count_g = count_t = 0
        # go to next nucleotide position
        nucleotide += 1

    return conserved_list
