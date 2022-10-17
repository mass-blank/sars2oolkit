# ALLELECOUNT SCRIPT
# Written by Nathaniel Lim (May 2014) for Joshua Mell and Rosie Redfield.
# Written for Python 3.4. Works with Python 3.3.1.  Does not work with Python 2.x.
# SNP counter script written for parsing samtools mpileup output.
# Intended to take samtools mpileup output as STDIN and gives count columns as STDOUT
# Use with samtools mpileup WITHOUT indicating a reference sequence (no -f option)
# Recommended usage: samtools mpileup -ABE -d 100000 -q1 reads.ref.bam | allelecount.v4.py | gzip > reads.ref.counts.txt.gz
# Upfront filtering using samtools view or sambamba view for finer control.
# Use of -q 1 option excludes multiply mapping reads, if using bwa mem alignments
# Output columns are: chr pos A C G T N a c g t n IN DEL indel
# Last field is only written for positions that indels were detected (i.e. leaves last column ragged).
# Indel field specifies counts of each insertion detected following the position and each deletion including the position.
# Revised on 2014-05-15. (rev 3)
# Revised again on 2016-03-14 (rev 4) to re-order ATCGN to ACGTN.

# Module Imports
import sys

# ----------------------------------------------------------------------------------------------------
# Base Counting Subroutine *[Completed]


def Base_Counter(InputRow):
    InputList = []
    InputList = InputRow.split(sep='\t')

    # Cleaning up Base String + Indel Counting
    CleanString = ''
    countIn = 0
    countDel = 0
    IndelHolder = []
    IndelDeterminant = 0
    CleanBool = False

    for currentIndex, Strholder in enumerate(InputList[4]):
        # Skipping of '^' Signage
        if CleanBool is True:
            CleanBool = False
            continue

        if Strholder == '^':
            CleanBool = True
            continue

        # Skipping Indel
        if IndelDeterminant > 0:
            IndelDeterminant -= 1
            continue

        if Strholder == '+':
            countIn += 1

            # Determining Length of Indel
            # Since Illumina NGS has an upper limit of less than 999bp read length
            IndelDeterminant = 0

            if (currentIndex + 4) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 4].isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 4]) + 3
                elif InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 2]) + 1
            elif (currentIndex + 3) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 2]) + 1

            IndelHolder.append(
                InputList[4][currentIndex:currentIndex + IndelDeterminant + 1])
            continue

        if Strholder == '-':
            countDel += 1

            # Determining Length of Indel
            # Since Illumina NGS has an upper limit of less than 999bp read length
            IndelDeterminant = 0

            if (currentIndex + 4) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 4].isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 4]) + 3
                elif InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 2]) + 1
            elif (currentIndex + 3) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() is True:
                    IndelDeterminant = int(
                        InputList[4][currentIndex + 1: currentIndex + 2]) + 1

            IndelHolder.append(
                InputList[4][currentIndex: currentIndex + len(str(IndelDeterminant)) + 1])
            continue

        CleanString += Strholder

    else:
        # Transferring Back Cleaned String
        InputList[4] = CleanString

        # '$' Signage Stripping
        InputList[4] = InputList[4].replace('$', '')

    # Base Count Var Initialization
    bigA = bigC = bigG = bigT = smallA = smallC = smallG = smallT = delBase = 0

    # Base Counting
    bigA = InputList[4].count('A')
    bigC = InputList[4].count('C')
    bigG = InputList[4].count('G')
    bigT = InputList[4].count('T')
    bigN = InputList[4].count('N')
    DOT = InputList[4].count('.')

    smallA = InputList[4].count('a')
    smallC = InputList[4].count('c')
    smallG = InputList[4].count('g')
    smallT = InputList[4].count('t')
    smallN = InputList[4].count('n')
    COMMA = InputList[4].count(',')

    delBase = InputList[4].count('*')

    # Internal Check - Throws Out Error (NOT STD-IN/OUT Compatible: Should break pipeline)
    InternalCounter = 0
    InternalCounter = bigA + bigC + bigG + bigT + bigN + smallA + \
        smallC + smallG + smallT + smallN + delBase + DOT + COMMA
    reported_number = int(InputList[3])
    if InternalCounter != reported_number and (not (reported_number == 0 and delBase == 1)):
        print('Error at position: ' + InputList[1])
        print('Reported count: ' + InputList[3])
        print('Internal counter: ' + str(InternalCounter))
        print('Internal sum: ' + str(len(InputList[4])))
        print('Number of insertions: ' + str(countIn))
        print('Number of deleions: ' + str(countDel))
        print('Post-processed bases: ' + InputList[4])
        sys.exit()

    # Indel Compilation
    IndelSetDict = set(IndelHolder)
    tmpIndelString = ''
    FinalIndelHolder = []

    for EveryIndel in IndelSetDict:
        tmpIndelString = ''
        tmpIndelString = str(IndelHolder.count(EveryIndel)) + ":" + EveryIndel
        FinalIndelHolder.append(tmpIndelString)

    # Return Output
    FinalOutput = InputList[0] + '	' + InputList[1] + '	' + str(bigA) + '	' + str(bigC) + '	' + str(bigG) + '	' + str(bigT) + '	' + str(bigN) + '	' + str(DOT) + '	' + str(
        smallA) + '	' + str(smallC) + '	' + str(smallG) + '	' + str(smallT) + '	' + str(smallN) + '	' + str(COMMA) + '	' + str(countIn) + '	' + str(countDel) + '	' + ';'.join(FinalIndelHolder)

    return FinalOutput


def read_pileup_write_allele(accession, input, output):
    try:
        with open(input, "r") as infile, open(output, "w") as outfile:
            for line in infile:
                generate_base = Base_Counter(line.strip())
                outfile.write(generate_base + "\n")
        print(f"{accession}: pileup created")
    except FileNotFoundError as ex:
        print(f"{ex}: File not found")

# Windows
# file_in = sys.argv[1]
# file_out = sys.argv[2]
# infile = open(file_in, "r")
# outfile = open(file_out, "w")
# for line in infile:
#     output = Base_Counter(line.strip())
#     outfile.write(output + "\n")
# infile.close()
# outfile.close()


# ----------------------------------------------------------------------------------------------------
# Running Caller
#StreamCollector = ''
# for EveryChar in sys.stdin.read():
#	if EveryChar == '\n':
#		print(Base_Counter(StreamCollector))
#		StreamCollector = ''
#	else:
#		StreamCollector += EveryChar
