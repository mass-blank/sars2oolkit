import os
import shutil
import argparse

from traceback import print_exc
from functions import *
from allelecount import *
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--infile',
        help='Input SRA accession .txt file',
        type=argparse.FileType('r'),
        required=True)
    parser.add_argument(
        '-d',
        '--download',
        action="store_true",
        help='Downloads accession list and outputs whether SARS2-CoV-2 positive')
    parser.add_argument(
        '-c',
        '--check',
        action="store_true",
        help="Check if the SRA's are Positive or Negative")
    parser.add_argument(
        '-b',
        '--bowtie',
        action='store_true',
        help='Aligns the reads to the reference genome')
    parser.add_argument(
        '-v',
        '--variants',
        action='store_true',
        help='Generates variants and logs them in: SRR#_mutations.txt')
    parser.add_argument(
        '-a',
        '--alleles',
        type=str,
        help='You must pass the nucleotide range.   Example:    --alleles 12345-54321')
    parser.add_argument(
        '-D',
        '--delete',
        action="store_true"
    )
    args = parser.parse_args()

# GLOBALS
my_mutations_text_file = Path(args.infile.name + "_mutations.txt")


if args.infile and args.alleles:
    nt_start, nt_stop = split_range(args.alleles)

    with open(args.infile.name) as f:
        lines = [line.rstrip() for line in f]
        for accession in lines:

            # INPUT
            my_sam_mpileup_file = Path(f"{accession}_pileup.txt")
            # OUTPUT
            my_alleles_text_file = Path(f"{accession}_{args.alleles}_NT.txt")

            if my_alleles_text_file.is_file():
                # PRINT FORMATTED OUTPUT
                print(f"\n{accession} \tNT\tA\tC\tG\tT")

                infile = open(my_alleles_text_file, 'r')
                lines = [line.rstrip() for line in infile]

                for line in lines:
                    print(line)
                infile.close()

            elif my_sam_mpileup_file.is_file():
                # THIS WRITES THE RANGE TO FILE
                print('SAM mpileup file exists')
                infile = open(my_sam_mpileup_file, "r")
                outfile = open(my_alleles_text_file, "w")

                for line in infile:
                    output = Base_Counter(line.strip())
                    outfile.write(output + '\n')
            else:
                # GENERATE MPILEUP
                gen_pileup(accession, nt_start, nt_stop)
                print('Pileup created')

if args.infile:
    file_variant = open(my_mutations_text_file.name, 'w+')

    with open(args.infile.name) as f:
        lines = [line.rstrip() for line in f]
        for accession in lines:

            my_sra_dir = Path(accession + "/")
            my_sra_file = Path(accession + "/" + accession + ".sra")
            my_fastq_file = Path(accession + ".fastq.gz")
            my_fastq_1_file = Path(accession + "_1.fastq.gz")
            my_fastq_2_file = Path(accession + "_2.fastq.gz")
            my_json_file = Path(accession + ".json")
            my_sam_file = Path(accession + ".sam")
            my_bcf_file = Path(accession + ".bcf")
            my_html_file = Path(accession + ".html")
            my_bam_file = Path(accession + ".bam")

# DOWNLOAD: files, check if positive for SARS-CoV-2

            if args.download:
                if (my_fastq_1_file.is_file() and my_fastq_2_file.is_file()) and (my_fastq_file.is_file() is False and my_json_file.is_file() is False):
                    if is_full() is True:
                        shutil.rmtree(my_sra_dir)
                    else:
                        pass
                    fastv_func(accession, my_fastq_1_file, my_fastq_2_file)
                elif my_fastq_file.is_file() and (my_fastq_1_file.is_file() is False and my_fastq_2_file.is_file() is False and my_json_file.is_file() is False):
                    fastv_func(accession)
                elif my_sra_file.is_file() and my_fastq_1_file.is_file() is False and my_fastq_2_file.is_file() is False and my_fastq_file.is_file() is False:
                    fastq_func(my_sra_file)
                elif my_sra_file.is_file() is False and my_fastq_1_file.is_file() is True and my_json_file.is_file() is False:
                    fastv_func(accession, my_fastq_1_file, my_fastq_2_file)
                elif my_json_file.is_file() is False and my_fastq_1_file.is_file() is False and my_json_file.is_file() is False:
                    fetch_func(accession + ".sra")
                    fastq_func(my_sra_file)
                    fastv_func(accession, my_fastq_1_file, my_fastq_2_file)
# BOWTIE
            elif args.bowtie:
                if my_fastq_file.is_file() and my_sam_file.is_file() is False:
                    # fastqc_func(accession)
                    bow_tie(accession)
                    sam_tools_view(accession)
                    sam_tools_sort(accession)
                    sam_tools_index(accession)
                    # gen_pileup(accession)
                elif my_fastq_1_file.is_file() and my_fastq_2_file.is_file() and my_sam_file.is_file() is False:
                    # fastqc_func(accession, my_fastq_1_file.name, my_fastq_2_file.name)
                    bow_tie(accession, my_fastq_1_file.name, my_fastq_2_file.name)
                    sam_tools_view(accession)
                    sam_tools_sort(accession)
                    sam_tools_index(accession)

                elif my_fastq_file.is_file() and my_fastq_2_file.is_file() and my_sam_file.is_file():
                    print("files exist")
# CALL VARIANTS
            elif args.variants:

                if my_bcf_file.is_file() is False:
                    call_mutations(accession)
                    file_variant.write(view_mutations(accession).stdout)
                else:
                    print("This has already been generated in file:" + my_mutations_text_file.name)
# DELETE
            elif args.delete:
                os.remove(my_sam_file)
                os.remove(my_bcf_file)
# CHECK IF SRA IS POSITIVE OR NEGATIVE
            elif args.check:
                if my_json_file.is_file():
                    try:
                        print('Removing SRA directory')
                        shutil.rmtree(my_sra_dir)
                    except Exception as ex:
                        print('The type of exception is: ', ex.__class__.__name__)
                        print_exc()
                    if check_positive(accession) == "NEGATIVE":
                        print(accession + " is Negative: deleting")
                        os.remove(my_json_file)
                        if my_fastq_file.is_file():
                            os.remove(my_fastq_file)
                        else:
                            os.remove(my_fastq_1_file)
                            os.remove(my_fastq_2_file)
                        os.remove(my_html_file)
                    elif check_positive(accession) == "POSITIVE":
                        print(accession + " is Positive")
                        os.remove(my_json_file)
                        os.remove(my_bam_file)
                        try:
                            os.remove(my_sam_file)
                        except Exception as ex:
                            print(ex.__class__.__name__)
                        if my_fastq_file.is_file():
                            os.remove(my_fastq_file)
                        elif my_fastq_1_file.is_file() and my_fastq_2_file.is_file():
                            os.remove(my_fastq_2_file)
                        else:
                            pass
                    else:
                        print('Neither Negative nor Positive')
                else:
                    print(my_json_file.name + " doesn't exist. Generate using by downloading sra and converting using --infile")
            elif args.alleles:
                pass

    file_variant.close()
else:
    print("No arguments")
