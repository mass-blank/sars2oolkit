import os
import shutil
import argparse

from functions import *
from allelecount import *
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--infile',
        help='Input SRA accession .txt file',
        type=argparse.FileType('r'))
    parser.add_argument(
        '-s',
        '--single',
        help='Input single SRA',
        type=str)
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

my_mutations_text_file = Path(args.infile.name + "_mutations.txt")

if args.infile and args.alleles:

    allele = args.alleles
    alleles = allele.split(',')
    print(alleles)

    with open(args.infile.name) as f:
        lines = [line.rstrip() for line in f]
        for accession in lines:

            # INPUT
            my_sam_mpileup_file = Path(f"{accession}_pileup.txt")
            # OUTPUT
            my_alleles_text_file = Path(f"{accession}_NT.txt")

            if my_alleles_text_file.is_file():
                # PRINT FORMATTED OUTPUT
                print(f"\n{accession} \tNT\tA\tC\tG\tT\tN\ta\tc\tg\tt\tn\tdel\tdot\tcomma")

                with open(my_alleles_text_file, 'r') as infile:
                    lines = [line.rstrip() for line in infile]
                    for line in lines:
                        split_line = line.split('\t')
                        for a in alleles:
                            if split_line[1].startswith(a):
                                print(line)

            elif my_sam_mpileup_file.is_file():
                # THIS WRITES THE RANGE TO FILE
                with open(my_sam_mpileup_file, "r") as infile, open(my_alleles_text_file, "w") as outfile:
                    for line in infile:
                        output = Base_Counter(line.strip())
                        outfile.write(output + '\n')
                print(f'{accession}: pileup created')
            elif my_sam_mpileup_file.is_file() is False and my_alleles_text_file.is_file() is False:
                # GENERATE MPILEUP
                gen_pileup(accession)
                with open(my_sam_mpileup_file, "r") as infile, open(my_alleles_text_file, "w") as outfile:
                    for line in infile:
                        output = Base_Counter(line.strip())
                        outfile.write(output + '\n')
                print(f'{accession}: pileup created')
            else:
                pass
if args.single:
    my_sra_dir = Path(args.single + "/")
    my_sra_file = Path(args.single + "/" + args.single + ".sra")
    my_fastq_file = Path(args.single + ".fastq.gz")
    my_fastq_1_file = Path(args.single + "_1.fastq.gz")
    my_fastq_2_file = Path(args.single + "_2.fastq.gz")
    my_json_file = Path(args.single + ".json")
    my_sam_file = Path(args.single + ".sam")
    my_bcf_file = Path(args.single + ".bcf")
    my_html_file = Path(args.single + ".html")
    my_bam_file = Path(args.single + ".bam")
    file_variant = open(args.single + '_mutations.txt', 'w+')
    if my_sra_file.is_file() is False:
        fetch_func(args.single + '.sra')
        fastq_func(my_sra_file)
    elif my_fastq_file.is_file():
        fastv_func(args.single)
        bow_tie(args.single)
        sam_tools_view(args.single)
        sam_tools_sort(args.single)
        sam_tools_index(args.single)
        call_mutations(args.single)
        file_variant.write(view_mutations(args.single).stdout)
    elif my_fastq_1_file.is_file() and my_fastq_2_file.is_file():
        fastv_func(args.single, my_fastq_1_file, my_fastq_2_file)
        bow_tie(args.single, my_fastq_1_file, my_fastq_2_file)
        sam_tools_view(args.single)
        sam_tools_sort(args.single)
        sam_tools_index(args.single)
        call_mutations(args.single)
        file_variant.write(view_mutations(args.single).stdout)
    else:
        print('Exists')


if args.infile:
    with open(my_mutations_text_file.name, 'w+') as file_variant, open(args.infile.name, 'r') as file_accession:
        lines = [line.rstrip() for line in file_accession]
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
            my_bam_file_sorted = Path(accession + ".sorted.bam")
            my_bam_file_index = Path(accession + ".sorted.bam.bai")

        # DOWNLOAD: files, check if positive for SARS-CoV-2
            if args.download:
                if (my_fastq_1_file.is_file()
                    and my_fastq_2_file.is_file()
                        and my_json_file.is_file() is False):
                    if is_full():
                        shutil.rmtree(my_sra_dir)
                    else:
                        pass
                    fastv_func(accession, my_fastq_1_file, my_fastq_2_file)
                    # print(1)
                elif (my_fastq_file.is_file()
                        and my_json_file.is_file() is False):
                    fastv_func(accession)
                    # print(2)
                elif (my_fastq_file.is_file()
                      and my_fastq_1_file.is_file()
                      and my_sra_file.is_file() is False):
                    fastq_exists(accession)
                    # print(3)
                elif (my_json_file.is_file()
                      and my_fastq_1_file.is_file()
                      and my_fastq_file.is_file()):
                    print('JSON and FASTQ files exist')
                elif (my_sra_file.is_file()
                      and my_fastq_1_file.is_file() is False
                      and my_fastq_file.is_file() is False):
                    fastq_func(my_sra_file)
                    # print(4)
                elif (my_sra_file.is_file()
                      and my_fastq_file.is_file()
                      and my_fastq_1_file.is_file() is False):
                    fastv_func(accession)
                elif (my_sra_file.is_file()
                      and my_fastq_1_file.is_file()
                      and my_fastq_file.is_file()
                      and my_sam_file.is_file()):
                    continue
                elif (my_fastq_file.is_file() is False
                      and my_fastq_1_file.is_file() is False
                      and my_sra_file.is_file() is False):
                    fastq_exists(accession)
                else:
                    # print(5)
                    pass

        # BOWTIE
            elif args.bowtie:
                if my_bam_file.is_file() is False:
                    bow_tie(accession, my_fastq_1_file, my_fastq_2_file, my_fastq_file)
                    sam_tools_view(accession)
                    sam_tools_sort(accession)
                    sam_tools_index(accession)
                elif (my_fastq_file.is_file() and my_sam_file.is_file()) and (my_fastq_1_file.is_file() and my_sam_file.is_file()):
                    print("FASTQ and .SAM files already exist. Proceed to next step.")
                else:
                    pass

        # CALL VARIANTS
            elif args.variants:
                if os.stat(my_mutations_text_file).st_size > 0:
                    with open(my_mutations_text_file, 'r') as mutations_file:
                        for line in mutations:
                            print(line)
                elif my_bcf_file.is_file() is False:
                    call_mutations(accession)
                elif my_bcf_file.is_file():
                    file_variant.write(view_mutations(accession).stdout)
                else:
                    pass

        # DELETE
            elif args.delete:
                try:
                    os.remove(my_sam_file)
                    os.remove(my_bcf_file)
                except FileNotFoundError as ex:
                    print(ex)
        # CHECK IF SRA IS POSITIVE OR NEGATIVE
            elif args.check:
                if my_json_file.is_file():
                    try:
                        shutil.rmtree(my_sra_dir)
                        if check_positive(my_json_file) == "NEGATIVE":
                            print(accession + " is Negative: deleting")
                            delete_accession(args.infile.name, accession)
                            os.remove(my_json_file)
                            os.remove(my_html_file)
                            os.remove(my_sam_file)
                            os.remove(my_bam_file)
                            os.remove(my_bcf_file)
                            os.remove(my_bam_file_sorted)
                            os.remove(my_bam_file_index)
                            if my_fastq_file.is_file():
                                os.remove(my_fastq_file)
                            else:
                                os.remove(my_fastq_1_file)
                                os.remove(my_fastq_2_file)
                        elif check_positive(my_json_file) == "POSITIVE" and mean_depth(my_json_file) <= 20.0:
                            print(f'Removing file with low depth < 15: {accession}')
                            delete_accession(args.infile.name, accession)
                            os.remove(my_json_file)
                            os.remove(my_html_file)
                            os.remove(my_sam_file)
                            os.remove(my_bam_file)
                            os.remove(my_bcf_file)
                            os.remove(my_bam_file_sorted)
                            os.remove(my_bam_file_index)
                        else:
                            pass
                    except FileNotFoundError as ex:
                        print(ex)
                else:
                    pass
