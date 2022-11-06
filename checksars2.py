import argparse
import os
import shutil
import xml.etree.ElementTree as ElTr
from collections import defaultdict
from pathlib import Path

import pandas as pd

from SRAClass.SRAClass import Accession
from SRAFunctions import functions
from SRAFunctions.conserved import conserved


class Bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def remove_files(args_file_name: str, accession: Accession):
    """This function removes all files associated with each accession,
    as defined in SRAClass.SRAClass"""
    functions.delete_accession(
        args_file_name, accession.acc)
    try:
        shutil.rmtree(accession.my_sra_dir)
    except OSError:
        print(f"{accession.my_sra_dir} doesn't exist")
    try:
        os.remove(accession.my_json_file)
    except FileNotFoundError:
        print(f"{accession.my_json_file} doesn't exist")
    try:
        os.remove(accession.my_html_file)
    except FileNotFoundError:
        print(f"{accession.my_html_file} doesn't exist")
    try:
        os.remove(accession.my_sam_file)
    except FileNotFoundError:
        print(f"{accession.my_sam_file} doesn't exist")
    try:
        os.remove(accession.my_bam_file)
    except FileNotFoundError:
        print(f"{accession.my_bam_file} doesn't exist")
    try:
        os.remove(accession.my_bcf_file)
    except FileNotFoundError:
        print(f"{accession.my_bcf_file} doesn't exist")
    try:
        os.remove(accession.my_bam_file_sorted)
    except FileNotFoundError:
        print(f"{accession.my_bam_file_sorted} doesn't exist")
    try:
        os.remove(accession.my_bam_file_index)
    except FileNotFoundError:
        print(f"{accession.my_bam_file_index} doesn't exist")
    if accession.my_fastq_file.exists():
        os.remove(accession.my_fastq_file)
    else:
        try:
            os.remove(accession.my_fastq_1_file)
        except FileNotFoundError:
            print(f"{accession.my_json_file} doesn't exist")
        try:
            os.remove(accession.my_fastq_2_file)
        except FileNotFoundError:
            print(f"{accession.my_fastq_2_file} doesn't exist")


def main():
    print(functions.progress())
    if args.infile.name is not None:
        my_mutations_text_file = Path(args.infile.name + "_mutations.txt")
    else:
        print('File not specified')

    if args.infile and args.alleles:
        # This block of code relates to the arguments: checksars2.py -i {file} -a
        # function from SRAFunctions/conserved.py returns a list of ancestral positions and corresponding allele
        alleles = dict(conserved())
        # A basal mutation that is excluded by conserved.py but is included here as it's been phylogenetically accepted
        # as an important mutation separating 2 lineages: v1 and a1
        # 3 = t as in [ 'a', 'c', 'g', 't' ]
        alleles[18060] = 3
        # initialize empty dictionary of list objects
        data = defaultdict(list)
        # creates a set of unique accessions
        my_acc_file_lines = set(functions.open_file_return_lines(args.infile.name))
        for idx, accession in enumerate(my_acc_file_lines):
            # creates an Accession object which creates Path objects
            acc = Accession(accession)
            print(f"{str(idx)}/{str(len(my_acc_file_lines))}", acc.acc)
            if acc.my_sam_mpileup_file.exists() is False and acc.my_bam_file_sorted.exists():
                functions.gen_pileup(acc)
                functions.read_pileup_write_allele(
                    acc.acc, acc.my_sam_mpileup_file, acc.my_alleles_text_file)
                my_allele_lines = functions.open_file_return_lines(acc.my_alleles_text_file)
                split_line = [line.split("\t") for line in my_allele_lines]
                row_dict = {int(item[1]): {'a': int(item[2]), 'c': int(item[3]), 'g': int(item[4]), 't': int(item[5])}
                            for item in split_line}
                for key in alleles.keys():
                    try:
                        percentages = functions.calculate_noise_return_percentages(row_dict[key]['a'],
                                                                                   row_dict[key]['c'],
                                                                                   row_dict[key]['g'],
                                                                                   row_dict[key]['t'])
                        data[acc.acc].append(percentages[alleles[key]])
                    except KeyError as ex:
                        print(f"Key doesn't exist {ex}: position {key}")

            elif acc.my_sam_mpileup_file.exists() and acc.my_alleles_text_file.exists() is False:
                # THIS WRITES THE RANGE TO FILE
                functions.read_pileup_write_allele(
                    acc.acc, acc.my_sam_mpileup_file, acc.my_alleles_text_file)
            elif acc.my_alleles_text_file.exists() and acc.my_sam_mpileup_file.exists():
                my_allele_lines = functions.open_file_return_lines(
                    acc.my_alleles_text_file)
                split_line = [line.split("\t") for line in my_allele_lines]
                row_dict = {int(item[1]): {'a': int(item[2]), 'c': int(item[3]), 'g': int(item[4]), 't': int(item[5])}
                            for item in split_line}
                for key in alleles.keys():
                    try:
                        percentages = functions.calculate_noise_return_percentages(row_dict[key]['a'],
                                                                                   row_dict[key]['c'],
                                                                                   row_dict[key]['g'],
                                                                                   row_dict[key]['t'])
                        data[acc.acc].append(percentages[alleles[key]])
                    except KeyError as ex:
                        print(f"Value {ex} doesn't exist at position {key}")
        df = pd.DataFrame.from_dict(data, orient="index")
        df = df.sum(axis=1)
        df = df.sort_values(axis=0)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            print(df)

    if args.infile and not args.alleles:
        lines = set(functions.open_file_return_lines(args.infile.name))
        for idx, accession in enumerate(lines):
            acc = Accession(accession)
            print(acc.acc)
            print(f"{idx + 1}/{len(lines)}")
            # DOWNLOAD: files, check if positive for SARS-CoV-2
            if args.download:

                if acc.my_fastq_1_file.exists() and acc.my_fastq_2_file.exists() and \
                        acc.my_json_file.exists() is False:
                    functions.fastv_func(acc.acc, acc.my_fastq_1_file, acc.my_fastq_2_file)

                elif acc.my_fastq_file.exists() and acc.my_json_file.exists() is False:
                    functions.fastv_func(acc.acc)

                elif acc.my_json_file.exists() and acc.my_fastq_1_file.exists() or \
                        acc.my_fastq_file.exists():
                    print("JSON and FASTQ files exist")

                elif acc.my_sra_file.exists() and (acc.my_fastq_file.exists() is False or
                                                   acc.my_fastq_1_file.exists() is False):
                    functions.fastq_func(acc)

                elif (acc.my_fastq_file.exists() is False or acc.my_fastq_1_file.exists() is False) and \
                        acc.my_bam_file_sorted.exists() is False:
                    print('fastq_exists function')

                    functions.fastq_exists(acc.acc)

                else:
                    pass

            # BOWTIE
            elif args.bowtie:
                if acc.my_bam_file.exists() is False and \
                        (acc.my_fastq_1_file.exists() or acc.my_fastq_file.exists()):
                    print(1)
                    functions.bow_tie(acc)
                    functions.sam_tools_view(acc)
                    functions.sam_tools_sort(acc)
                    functions.sam_tools_index(acc)
                elif acc.my_bam_file.exists() is False and \
                        (acc.my_fastq_file.exists() or acc.my_fastq_1_file.exists()):
                    functions.sam_tools_view(acc)
                    functions.sam_tools_sort(acc)
                    functions.sam_tools_index(acc)
                elif acc.my_bam_file.exists() and acc.my_bam_file_sorted.exists():
                    print("Alignment complete")
                elif acc.my_fastq_1_file.exists() is False or acc.my_fastq_file.exists() is False:
                    print("No .FASTQ files")
                else:
                    pass

            # CALL VARIANTS
            elif args.variants:
                if os.stat(my_mutations_text_file).st_size > 0:
                    with open(my_mutations_text_file, "r") as mutations_file:
                        for line in mutations_file:
                            print(line)
                elif acc.my_bcf_file.exists() is False:
                    functions.call_mutations(acc)
                elif acc.my_bcf_file.exists():
                    with open(my_mutations_text_file, 'w+') as file_variant:
                        file_variant.write(
                            functions.view_mutations(acc).stdout)
                else:
                    pass

            # DELETE OXFORD NANOPORE
            elif args.delete:
                try:
                    result = functions.xml_parse(acc.acc, 'Platform')[0]
                except ElTr.ParseError as ex:
                    print(ex)
                finally:
                    if result == 'OXFORD_NANOPORE':
                        remove_files(args.infile.name, acc)

            # CHECK IF SRA IS POSITIVE OR NEGATIVE
            elif args.check:
                if acc.my_json_file.exists():
                    depth = functions.mean_depth(acc.my_json_file)
                    is_positive = functions.is_positive(acc.my_json_file)
                    if not is_positive:
                        print(acc.acc + " is Negative: deleting")
                        remove_files(args.infile.name, acc)
                    elif is_positive and (depth <= args.depth):
                        remove_files(args.infile.name, acc)
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
                        help="Downloads accession list and if accession is SARS-CoV-2 positive and depth"
                        )
    parser.add_argument("-c", "--check",
                        action="store_true",
                        help="Check if the SRA's are Positive or Negative and delete accession/files"
                        )
    parser.add_argument("-T", "--depth",
                        type=float, default=100.00,
                        help="Mean depth for acceptable positive"
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
