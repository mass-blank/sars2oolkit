import argparse
import os
import shutil
from collections import defaultdict
from pathlib import Path
from sre_constants import RANGE

import numpy as np
import pandas as pd

from allelecount import *
from conserved import conserved
from functions import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--infile",
        help="Input SRA accession .txt file",
        type=argparse.FileType("r"),
    )
    parser.add_argument(
        "-d",
        "--download",
        action="store_true",
        help="Downloads accession list and outputs whether SARS2-CoV-2 positive",
    )
    parser.add_argument(
        "-c",
        "--check",
        action="store_true",
        help="Check if the SRA's are Positive or Negative",
    )
    parser.add_argument(
        "--depth", type=float, default=15.00, help="Mean depth for acceptable positive."
    )
    parser.add_argument(
        "-b",
        "--bowtie",
        action="store_true",
        help="Aligns the reads to the reference genome",
    )
    parser.add_argument(
        "-v",
        "--variants",
        action="store_true",
        help="Generates variants and logs them in: SRR#_mutations.txt",
    )
    parser.add_argument(
        "-a",
        "--alleles",
        action="store_true",
        # type=str,
        help="You must pass the nucleotides seperated by commas.   Example:    --alleles 12345,54321",
    )
    parser.add_argument(
        "-t", "--total", action="store_true", help="generates all alleles"
    )
    parser.add_argument("-D", "--delete", action="store_true")

    args = parser.parse_args()


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


class Accession:
    def __init__(self, accession):
        self.my_sra_dir = Path(f"{accession}/")
        self.my_sra_file = Path(f"{accession}/{accession}.sra")
        self.my_fastq_file = Path(accession + ".fastq.gz")
        self.my_fastq_1_file = Path(accession + "_1.fastq.gz")
        self.my_fastq_2_file = Path(accession + "_2.fastq.gz")
        self.my_json_file = Path(accession + ".json")
        self.my_sam_file = Path(accession + ".sam")
        self.my_bcf_file = Path(accession + ".bcf")
        self.my_html_file = Path(accession + ".html")
        self.my_bam_file = Path(accession + ".bam")
        self.my_bam_file_sorted = Path(accession + ".sorted.bam")
        self.my_bam_file_index = Path(accession + ".sorted.bam.bai")
        self.my_alleles_text_file = Path(f"{accession}_NT.txt")
        self.my_sam_mpileup_file = Path(f"{accession}_pileup.txt")


ACC_RANGE = 1085

if args.infile.name is not None:
    my_mutations_text_file = Path(args.infile.name + "_mutations.txt")

if args.total:
    with open(args.infile.name) as f:
        lines = [line.rstrip() for line in f]
        for accession in lines:
            acc = Accession(accession)
            if acc.my_alleles_text_file.is_file():
                # PRINT FORMATTED OUTPUT
                print(
                    f"\n{accession} \tNT\tA\tC\tG\tT\tN\ta\tc\tg\tt\tn\tdel\tdot\tcomma"
                )

                with open(acc.my_alleles_text_file, "r") as infile:
                    lines = [line.rstrip() for line in infile]
                    for line in lines:
                        print(line)

if args.infile and args.alleles:

    # alleles = args.alleles.split(',')
    alleles = conserved()
    # print(alleles)
    legend = ["a", "c", "g", "t"]
    data = defaultdict(list)

    with open(args.infile.name) as infile:
        read_lines = [line_infile.rstrip() for line_infile in infile]
        for idx, accession in enumerate(read_lines[0:ACC_RANGE]):
            acc = Accession(accession)
            if (
                acc.my_sam_mpileup_file.is_file() is False
                and acc.my_alleles_text_file.is_file() is False
            ):
                # GENERATE MPILEUP
                gen_pileup(accession)
                read_pileup_write_allele(
                    accession, acc.my_sam_mpileup_file, acc.my_alleles_text_file
                )
            elif acc.my_sam_mpileup_file.is_file():
                # THIS WRITES THE RANGE TO FILE
                read_pileup_write_allele(
                    accession, acc.my_sam_mpileup_file, acc.my_alleles_text_file
                )
            elif acc.my_alleles_text_file.is_file():
                # PRINT FORMATTED OUTPUT
                print(f"{str(idx)}/{str(len(read_lines[0:ACC_RANGE]))}")
                # print(
                #     f"\n{accession} \tNT\tA\tC\tG\tT\tN\ta\tc\tg\tt\tn\tdel\tdot\tcomma"
                # )
                try:
                    with open(acc.my_alleles_text_file, "r") as a_file:
                        lines = [line.rstrip() for line in a_file]
                        for line in lines:
                            split_line = line.split("\t")
                            for allele in alleles:
                                if int(split_line[1]) == allele[0]:
                                    nt_array = np.array(
                                        [
                                            int(split_line[2]),
                                            int(split_line[3]),
                                            int(split_line[4]),
                                            int(split_line[5]),
                                        ]
                                    )
                                    percentages = nt_array / nt_array.sum(axis=0)
                                    # print(line)
                                    # print(
                                    # f"\t\t{bcolors.OKGREEN}{allele[0]}{str.capitalize(legend[allele[1]])}{bcolors.ENDC}\t{percentages[0]:.2%}\t{percentages[1]:.2%}\t{percentages[2]:.2%}\t{percentages[3]:.2%}\n"
                                    # )
                                    data[accession].append(percentages[allele[1]])
                except FileNotFoundError as ex:
                    print(f"{ex} File not found")
            else:
                continue

    df = pd.DataFrame.from_dict(data, orient="index")
    df = df.sum(axis=1)
    df = df.sort_values(axis=0)
    print(df)


if args.infile:
    with open(my_mutations_text_file.name, "w+") as file_variant, open(
        args.infile.name, "r"
    ) as file_accession:
        lines = [line.rstrip() for line in file_accession]
        for accession in lines:
            acc = Accession(accession)
            # DOWNLOAD: files, check if positive for SARS-CoV-2
            if args.download:
                if (
                    acc.my_fastq_1_file.is_file()
                    and acc.my_fastq_2_file.is_file()
                    and acc.my_json_file.is_file() is False
                ):
                    if is_full():
                        shutil.rmtree(acc.my_sra_dir)
                    else:
                        pass
                    fastv_func(accession, acc.my_fastq_1_file, acc.my_fastq_2_file)
                elif (
                    acc.my_fastq_file.is_file() and acc.my_json_file.is_file() is False
                ):
                    # fastv_func(accession)
                    continue
                elif (
                    acc.my_fastq_file.is_file()
                    and acc.my_fastq_1_file.is_file()
                    and acc.my_sra_file.is_file() is False
                ):
                    fastq_exists(accession)
                elif (
                    acc.my_json_file.is_file()
                    and acc.my_fastq_1_file.is_file()
                    and acc.my_fastq_file.is_file()
                ):
                    print("JSON and FASTQ files exist")
                elif (
                    acc.my_sra_file.is_file()
                    and acc.my_fastq_1_file.is_file() is False
                    and acc.my_fastq_file.is_file() is False
                ):
                    fastq_func(acc.my_sra_file)
                elif (
                    acc.my_sra_file.is_file()
                    and acc.my_fastq_file.is_file()
                    and acc.my_fastq_1_file.is_file() is False
                ):
                    # fastv_func(accession)
                    continue
                elif (
                    acc.my_sra_file.is_file()
                    and acc.my_fastq_1_file.is_file()
                    and acc.my_fastq_file.is_file()
                    and acc.my_sam_file.is_file()
                ):
                    continue
                elif (
                    acc.my_fastq_file.is_file() is False
                    and acc.my_fastq_1_file.is_file() is False
                    and acc.my_sra_file.is_file() is False
                ):
                    fastq_exists(accession)
                else:
                    pass

            # BOWTIE
            elif args.bowtie:
                if acc.my_bam_file.is_file() is False:
                    bow_tie(
                        accession,
                        acc.my_fastq_1_file,
                        acc.my_fastq_2_file,
                        acc.my_fastq_file,
                    )
                    sam_tools_view(accession)
                    sam_tools_sort(accession)
                    sam_tools_index(accession)
                elif (acc.my_fastq_file.is_file() and acc.my_sam_file.is_file()) and (
                    acc.my_fastq_1_file.is_file() and acc.my_sam_file.is_file()
                ):
                    print("FASTQ and .SAM files already exist. Proceed to next step.")
                else:
                    pass

            # CALL VARIANTS
            elif args.variants:
                if os.stat(my_mutations_text_file).st_size > 0:
                    with open(my_mutations_text_file, "r") as mutations_file:
                        for line in mutations_file:
                            print(line)
                elif acc.my_bcf_file.is_file() is False:
                    call_mutations(accession)
                elif acc.my_bcf_file.is_file():
                    file_variant.write(view_mutations(accession).stdout)
                else:
                    pass

            # DELETE
            elif args.delete:
                try:
                    os.remove(acc.my_sam_file)
                    os.remove(acc.my_bcf_file)
                except FileNotFoundError as ex:
                    print(ex)
            # CHECK IF SRA IS POSITIVE OR NEGATIVE
            elif args.check:
                if acc.my_json_file.is_file():
                    try:
                        shutil.rmtree(acc.my_sra_dir)
                        if check_positive(acc.my_json_file) == "NEGATIVE":
                            print(accession + " is Negative: deleting")
                            delete_accession(args.infile.name, accession)
                            os.remove(acc.my_json_file)
                            os.remove(acc.my_html_file)
                            os.remove(acc.my_sam_file)
                            os.remove(acc.my_bam_file)
                            os.remove(acc.my_bcf_file)
                            os.remove(acc.my_bam_file_sorted)
                            os.remove(acc.my_bam_file_index)
                            if acc.my_fastq_file.is_file():
                                os.remove(acc.my_fastq_file)
                            else:
                                os.remove(acc.my_fastq_1_file)
                                os.remove(acc.my_fastq_2_file)
                        elif (
                            check_positive(acc.my_json_file) == "POSITIVE"
                            and mean_depth(acc.my_json_file) <= args.depth
                        ):
                            print(
                                f"Removing file with low depth < {mean_depth}: {accession}"
                            )
                            delete_accession(args.infile.name, accession)
                            os.remove(acc.my_json_file)
                            os.remove(acc.my_html_file)
                            os.remove(acc.my_sam_file)
                            os.remove(acc.my_bam_file)
                            os.remove(acc.my_bcf_file)
                            os.remove(acc.my_bam_file_sorted)
                            os.remove(acc.my_bam_file_index)
                        else:
                            pass
                    except FileNotFoundError as ex:
                        print(ex)
                else:
                    pass
