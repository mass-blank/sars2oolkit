import argparse
import os
import shutil
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import functions
from allelecount import read_pileup_write_allele
from conserved import conserved


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
        self.acc = accession
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


class Conserved:
    def __init__(self, position_allele):
        self.nucleotide: int = position_allele[0]
        self.allele = position_allele[1]

    def get_allele(self, pos):
        self.legend = ["a", "c", "g", "t"]
        return self.legend[pos]


class Line:
    def __init__(self, nucleotide, bigA, bigC, bigG, bigT):
        self.nucleotide = int(nucleotide)
        self.A: int = int(bigA)
        self.C: int = int(bigC)
        self.G: int = int(bigG)
        self.T: int = int(bigT)


def calculate_noise_return_percentages(rowA, rowC, rowG, rowT):
    nt_array = np.array(
        [rowA, rowC, rowG, rowT])
    A, B = np.partition(nt_array, 1)[0:2]
    noise = A + B / 2
    nt_array = nt_array - noise
    nt_array = nt_array.clip(min=0)
    perc_array = nt_array / nt_array.sum(axis=0)
    return perc_array


def open_file_return_lines(file) -> list:
    """Open file return lines as a list stripped of \\n lines"""
    try:
        with open(file, 'r') as f:
            lines = [line.rstrip() for line in f]
        return lines
    except FileNotFoundError as ex:
        print(f'File not found: {ex}')


def main():
    ACC_RANGE = 1618

    if args.infile is not None:
        my_mutations_text_file = Path(args.infile.name + "_mutations.txt")
    else:
        print('File not specified')

    if args.infile and args.alleles:
        alleles = dict(conserved())
        data = defaultdict(list)
        my_acc_file_lines = set(open_file_return_lines(args.infile.name))
        for idx, accession in enumerate(my_acc_file_lines):
            acc = Accession(accession)
            print(f"{str(idx)}/{str(len(my_acc_file_lines))}", acc.acc)
            if acc.my_bam_file_sorted.is_file() is True and acc.my_sam_mpileup_file.is_file() is False:
                # GENERATE MPILEUP
                print(1)
                functions.gen_pileup(acc.acc)
                read_pileup_write_allele(
                    acc.acc, acc.my_sam_mpileup_file, acc.my_alleles_text_file)
                my_allele_lines = open_file_return_lines(
                    acc.my_alleles_text_file)
                for line in my_allele_lines:
                    split_line = line.split("\t")
                    allele_acc_row = Line(
                        split_line[1], split_line[2], split_line[3], split_line[4], split_line[5])
                    if allele_acc_row.nucleotide in alleles.keys():
                        percentages = calculate_noise_return_percentages(
                            allele_acc_row.A, allele_acc_row.C, allele_acc_row.G, allele_acc_row.T)
                        # print(line)
                        # print(
                        #     f"\t\t{bcolors.OKGREEN}{allele_acc_row.nucleotide}{str.capitalize(legend[alleles[allele_acc_row.nucleotide]])}{bcolors.ENDC}\t{percentages[0]:.2%}\t{percentages[1]:.2%}\t{percentages[2]:.2%}\t{percentages[3]:.2%}\n"
                        # )
                        # dictionary of accessions and noise per conserved nucleotide
                        data[acc.acc].append(
                            percentages[alleles[allele_acc_row.nucleotide]])
            elif acc.my_bam_file_sorted.is_file() and acc.my_sam_mpileup_file.is_file() and acc.my_alleles_text_file.is_file() is False:
                # THIS WRITES THE RANGE TO FILE
                print(2)
                read_pileup_write_allele(
                    acc.acc, acc.my_sam_mpileup_file, acc.my_alleles_text_file)
            elif acc.my_alleles_text_file.is_file():
                print(3)
                # print(
                #     f"\n{accession} \tNT\tA\tC\tG\tT\tN\ta\tc\tg\tt\tn\tdel\tdot\tcomma"
                # )
                my_allele_lines = open_file_return_lines(
                    acc.my_alleles_text_file)
                for line in my_allele_lines:
                    split_line = line.split("\t")
                    allele_acc_row = Line(
                        split_line[1], split_line[2], split_line[3], split_line[4], split_line[5])
                    if allele_acc_row.nucleotide in alleles.keys():
                        percentages = calculate_noise_return_percentages(
                            allele_acc_row.A, allele_acc_row.C, allele_acc_row.G, allele_acc_row.T)
                        # print(line)
                        # print(
                        #     f"\t\t{bcolors.OKGREEN}{allele_acc_row.nucleotide}{str.capitalize(legend[alleles[allele_acc_row.nucleotide]])}{bcolors.ENDC}\t{percentages[0]:.2%}\t{percentages[1]:.2%}\t{percentages[2]:.2%}\t{percentages[3]:.2%}\n"
                        # )
                        # dictionary of accessions and noise per conserved nucleotide
                        data[acc.acc].append(
                            percentages[alleles[allele_acc_row.nucleotide]])

        df = pd.DataFrame.from_dict(data, orient="index")
        df = df.sum(axis=1)
        df = df.sort_values(axis=0)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(df)

    if args.infile and not args.alleles:
        with open(my_mutations_text_file.name, "w+") as file_variant:
            lines = set(open_file_return_lines(args.infile.name))
            for idx, accession in enumerate(lines):
                acc = Accession(accession)
                print(f"{idx}/{len(lines)}")
                # DOWNLOAD: files, check if positive for SARS-CoV-2
                if args.download:
                    if (acc.my_fastq_1_file.is_file() and acc.my_fastq_2_file.is_file() and acc.my_json_file.is_file() is False):
                        if functions.is_full():
                            shutil.rmtree(acc.my_sra_dir)
                        else:
                            pass
                        # fastv_func(acc.acc, acc.my_fastq_1_file,
                            #    acc.my_fastq_2_file)
                    elif (acc.my_fastq_file.is_file() and acc.my_json_file.is_file() is False):
                        # fastv_func(acc.acc)
                        continue
                    elif (acc.my_fastq_file.is_file() and acc.my_fastq_1_file.is_file() and acc.my_sra_file.is_file() is False):
                        functions.fastq_exists(acc.acc)
                    elif (acc.my_json_file.is_file() and acc.my_fastq_1_file.is_file() and acc.my_fastq_file.is_file()):
                        print("JSON and FASTQ files exist")
                    elif (acc.my_sra_file.is_file() and acc.my_fastq_1_file.is_file() is False and acc.my_fastq_file.is_file() is False):
                        functions.fastq_func(acc.my_sra_file)
                    elif (acc.my_sra_file.is_file() and acc.my_fastq_file.is_file() and acc.my_fastq_1_file.is_file() is False):
                        # fastv_func(acc.acc)
                        continue
                    elif (acc.my_sra_file.is_file() and acc.my_fastq_1_file.is_file() and acc.my_fastq_file.is_file() and acc.my_sam_file.is_file()):
                        continue
                    elif (
                        acc.my_fastq_file.is_file() is False
                        and acc.my_fastq_1_file.is_file() is False
                        and acc.my_sra_file.is_file() is False
                    ):
                        functions.fastq_exists(acc.acc)
                    else:
                        pass

                # BOWTIE
                elif args.bowtie:
                    if acc.my_bam_file.is_file() is False:
                        functions.bow_tie(
                            acc.acc,
                            acc.my_fastq_1_file,
                            acc.my_fastq_2_file,
                            acc.my_fastq_file,
                        )
                        functions.sam_tools_view(acc.acc)
                        functions.sam_tools_sort(acc.acc)
                        functions.sam_tools_index(acc.acc)
                    elif (acc.my_fastq_file.is_file() and acc.my_sam_file.is_file()) or (
                        acc.my_fastq_1_file.is_file() and acc.my_sam_file.is_file()
                    ):
                        print(
                            "FASTQ and .SAM files already exist. Proceed to next step.")
                    else:
                        pass

                # CALL VARIANTS
                elif args.variants:
                    if os.stat(my_mutations_text_file).st_size > 0:
                        with open(my_mutations_text_file, "r") as mutations_file:
                            for line in mutations_file:
                                print(line)
                    elif acc.my_bcf_file.is_file() is False:
                        functions.call_mutations(acc.acc)
                    elif acc.my_bcf_file.is_file():
                        file_variant.write(
                            functions.view_mutations(acc.acc).stdout)
                    else:
                        pass

                # DELETE OXFORD NANOPORE
                elif args.delete:
                    try:
                        result = functions.xml_parse(acc.acc, 'Platform')[0]
                        if result == 'OXFORD_NANOPORE':
                            print(result)
                            functions.delete_accession(
                                args.infile.name, acc.acc)
                            try:
                                os.remove(acc.my_sra_file)
                                os.remove(acc.my_fastq_file)
                                os.remove(acc.my_json_file)
                                os.remove(acc.my_html_file)
                                os.remove(acc.my_sam_file)
                                os.remove(acc.my_bam_file)
                                os.remove(acc.my_bcf_file)
                                os.remove(acc.my_bam_file_sorted)
                                os.remove(acc.my_bam_file_index)
                            except FileNotFoundError as ex:
                                print(f"{ex}")
                    except ET.ParseError as ex:
                        print(ex)

                # CHECK IF SRA IS POSITIVE OR NEGATIVE
                elif args.check:
                    if acc.my_json_file.is_file():
                        try:
                            shutil.rmtree(acc.my_sra_dir)
                            if functions.check_positive(acc.my_json_file) == "NEGATIVE":
                                print(acc.acc + " is Negative: deleting")
                                functions.delete_accession(
                                    args.infile.name, acc.acc)
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
                                functions.check_positive(
                                    acc.my_json_file) == "POSITIVE"
                                and functions.mean_depth(acc.my_json_file) <= args.depth
                            ):
                                print(
                                    f"Removing file with low depth < {functions.mean_depth}: {acc.acc}")
                                functions.delete_accession(
                                    args.infile.name, acc.acc)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile",
                        help="Input SRA accession .txt file",
                        type=argparse.FileType("r")
                        )
    parser.add_argument("-d", "--download",
                        action="store_true",
                        help="Downloads accession list and outputs whether SARS2-CoV-2 positive"
                        )
    parser.add_argument("-c", "--check",
                        action="store_true",
                        help="Check if the SRA's are Positive or Negative"
                        )
    parser.add_argument("-T", "--depth",
                        type=float, default=20.00,
                        help="Mean depth for acceptable positive."
                        )
    parser.add_argument("-b", "--bowtie",
                        action="store_true",
                        help="Aligns the reads to the reference genome"
                        )
    parser.add_argument("-v", "--variants",
                        action="store_true",
                        help="Generates variants and logs them in: SRR#_mutations.txt"
                        )
    parser.add_argument("-a", "--alleles",
                        action="store_true",
                        help="This function generates alleles for conserved positions"
                        )
    parser.add_argument("-D", "--delete",
                        action="store_true"
                        )
    args = parser.parse_args()
    main()
