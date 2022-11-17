import json
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ElTr
from pathlib import Path
from typing import Generator

import numpy as np

from SRAClass.SRAClass import Accession
from SRAFunctions.allelecount import Base_Counter





def progress():
    """How many sorted.bam files exist out of total accession"""
    subdir = Path(__file__).parent.resolve()
    SRR_FILE = Path(subdir,'SRR_List_1.txt')
    lines = open_file_return_lines(SRR_FILE)
    count = 0
    for line in set(lines):
        accession = Accession(line)
        if accession.my_bam_file_sorted.exists():
            count += 1
    return count, len(lines)


def calculate_noise_return_percentages(row_a: int, row_c: int, row_g: int, row_t: int):
    """Calculate the noise in each row and return percentages(float) without noise"""

    nt_array = np.array(
        [row_a, row_c, row_g, row_t])
    if nt_array.sum() > 30:
        a, b = np.partition(nt_array, 1)[0:2]
        noise = a + b / 2
        nt_array = nt_array - noise
        nt_array = nt_array.clip(min=0)
        perc_array = nt_array / nt_array.sum(axis=0)
        return perc_array
    else:
        nt_array = np.array([0, 0, 0, 0])
        return nt_array


def open_file_return_lines(file: Path) -> list:
    with open(file, 'r') as f:
        lines = [line.rstrip() for line in f]
    return lines

def open_file_return_generator(file: Path) -> Generator:
    with open(file, 'r') as f:
        for line in f:
            yield line.rstrip()

def xml_parse(accession: str, search: str) -> list:
    """Returns bool for specific search item found in the xml file for an accession from NCBI"""
    args = f"esummary -db sra -id {accession}"
    myxml = subprocess.run(args, shell=True, capture_output=True, text=True)
    myroot = ElTr.fromstring(myxml.stdout)
    result = [result.text for result in myroot.findall(f'.//{search}')]
    return result


def call_mutations(accession: Accession):
    """calls mutations and outputs them into a bcf file"""
    args = f"bcftools mpileup -d 150000 -Ou -f SARS-CoV-2_reference.fasta {accession.my_bam_file_sorted} | " \
           f"bcftools call -mv -Ob -o {accession.my_bcf_file}"
    subprocess.run(args, shell=True)


def gen_pileup(accession: Accession):
    """Generate pileup"""
    args = (
        f"samtools mpileup -d 10000 {accession.my_bam_file_sorted} -o {accession.my_sam_mpileup_file} "
    )
    try:
        subprocess.run(args, shell=True)
    except subprocess.CalledProcessError as ex:
        print(f"{ex.output} Couldn't run mpileup")


def view_mutations(accession: Accession):
    """View the called mutations"""
    args = "bcftools query -f '" + accession.acc + \
           r" %REF%POS%ALT\n' " + accession.my_bcf_file
    return subprocess.run(args, shell=True, capture_output=True, text=True)


def dir_is_empty(path: Path) -> bool:
    """Is the directory empty returns bool"""
    if os.path.exists(path) and not os.path.isfile(path):
        if not os.listdir(path):
            return True
        else:
            return False
    else:
        pass


def glob_re(pattern, strings):
    return filter(re.compile(pattern).match, strings)


def one_or_two_fastq_gz(accession, re_result):
    """Checks if fastq file is Paired end or Single stranded"""
    current_dir = os.getcwd()
    if len(re_result) == 1:
        print(1)
        re_result = list(re_result)
        infile = Path(accession + "/" + re_result[0])
        outfile = Path(accession + ".fastq.gz")
        os.rename(infile, os.path.join(current_dir, outfile))
        print(f"{outfile.name} has finished downloading")
    elif len(re_result) == 2:
        print(2)
        print(re_result)
        r = list(sorted(re_result))
        file_1 = r[0]
        file_2 = r[1]
        print('File 1: ', file_1)
        print('File 2: ', file_2)
        infile1 = Path(accession + "/" + file_1)
        infile2 = Path(accession + "/" + file_2)
        outfile1 = Path(accession + "_1.fastq.gz")
        outfile2 = Path(accession + "_2.fastq.gz")
        os.rename(infile1, os.path.join(current_dir, outfile1))
        os.rename(infile2, os.path.join(current_dir, outfile2))
        print(f"{outfile1.name} and {outfile2.name} have finished downloading")


def fastq_exists(accession):
    """Checks if fastq files already exist in the SRA folder if not downloads them"""
    if dir_is_empty(accession) is False:
        filenames = glob_re(r'.*\.f.*(gz)?', os.listdir(accession))
        results = [item for item in filenames]
        results.sort()
        one_or_two_fastq_gz(accession, results)
    else:
        result = subprocess.run(
            ["prefetch", "--type", "fastq", accession], capture_output=True
        )
        if result.stderr.decode("utf-8").count("\n") == 2:
            subprocess.run(["prefetch", accession])
        else:
            re_result = re.findall(r"'.*\.f.*[gz]?'", result.stderr.decode("utf-8"))
            re_result = [i.strip("'") for i in re_result]
            re_result = set(re_result)
            one_or_two_fastq_gz(accession, re_result)


def sra_is_paired(sra_file):
    if sra_file.exists():
        contents = subprocess.check_output(
            ["fastq-dump", "-X", "1", "-Z", "--split-spot", sra_file]
        )
        if contents.count(b"\n") == 4:
            return False
        elif contents.count(b"\n") == 8:
            return True
        else:
            pass
    else:
        print(f"{sra_file} does not exist.")


def bow_tie(accession: Accession):
    """Aligns fastq files with bowtie2 if paired end reads or bwa if single stranded"""
    if accession.my_fastq_1_file.exists() and accession.my_fastq_2_file.exists():
        args = f"bowtie2 -p 8 -x bowtie -1 {accession.my_fastq_1_file} -2 {accession.my_fastq_2_file} " \
               f"-S {accession.my_sam_file}"
        subprocess.run(args, shell=True)
    elif accession.my_fastq_file.exists():
        args = f"bwa mem -t 8 SARS-CoV-2_reference.fasta {accession.my_fastq_file} > {accession.my_sam_file}"
        subprocess.run(args, shell=True)
    else:
        print("No FASTQ files.")


def sam_tools_view(accession: Accession):
    """Converts sam file to bam file"""
    args = f"samtools view -@ 8 -bS {accession.my_sam_file} > {accession.my_bam_file}"
    subprocess.run(args, shell=True)


def sam_tools_sort(accession: Accession):
    """Converts bam file to sorted.bam file"""
    args = f"samtools sort -@ 8 {accession.my_bam_file} -o {accession.my_bam_file_sorted}"
    subprocess.run(args, shell=True)


def sam_tools_index(accession: Accession):
    """Creates an index file for sorted.bam file"""
    args = f"samtools index -@ 8 {accession.my_bam_file_sorted}"
    subprocess.run(args, shell=True)


def fastq_func(accession: Accession):
    """Downloads either split fastq files or single stranded"""
    if sra_is_paired(accession.my_sra_file):
        args = f"parallel-fastq-dump --sra-id {accession.my_sra_file} --threads 4 --gzip"
        # args = f"fastq-dump --gzip {accession.my_sra_file}"
        subprocess.run(args, shell=True)
    else:
        # args = f"parallel-fastq-dump --sra-id {accession.my_sra_file} --threads 4 --split-files --gzip"
        args = f"fastq-dump --split-3 --gzip {accession.my_sra_file}"
        subprocess.run(args, shell=True)


def fetch_func(accession):
    subprocess.run(["prefetch", accession], stdout=subprocess.PIPE)


def fastqc_func(accession, fastq_file_1=None, fastq_file_2=None):
    if sra_is_paired(accession):
        args = f"fastqc {fastq_file_1} {fastq_file_2}"
        subprocess.run(args, shell=True)
    else:
        args = f"fastqc {accession}.fastq.gz"
        subprocess.run(args, shell=True)


def fastv_func(accession, fastq_file_1=None, fastq_file_2=None):
    """Creates a fastv report"""
    if fastq_file_1 is not None and fastq_file_2 is not None:
        if fastq_file_1.exists() and fastq_file_2.exists():
            args = f" fastv --in1 {fastq_file_1} --in2 {fastq_file_2} -g SARS-CoV-2.genomes.fa " \
                   f"-k SARS-CoV-2.kmer.fa -h {accession}.html -j {accession}.json"
            subprocess.run(args, shell=True)
    else:
        args = f"fastv -i {accession}.fastq.gz -g SARS-CoV-2.genomes.fa -k SARS-CoV-2.kmer.fa " \
               f"-h {accession}.html -j {accession}.json"
        subprocess.run(args, shell=True)


def is_positive(json_file) -> bool:
    """Checks if an accession is positive for SARS-CoV-2 (JSON file)"""
    with open(json_file, "r") as file:
        data = json.loads(file.read())
        if data["kmer_detection_result"]["result"] == "POSITIVE":
            return True
        else:
            return False


def delete_accession(file, accession):
    """Deletes accession from accession file"""
    with open(file, "r+") as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            if accession not in line.strip("\n"):
                file.write(line)
        file.truncate()


def mean_depth(json_file):
    """Returns the mean coverage of an accession"""
    with open(json_file, "r") as file:
        data = json.loads(file.read())
        return float(data["kmer_detection_result"]["mean_coverage"])


def is_full():
    total, used, free = shutil.disk_usage("/")
    free = free / 2 ** 20
    if free < 900:
        print(f"There is only {free} free on the hard disk. Deleting...")
        return True
    else:
        return False


def read_pileup_write_allele(accession, input_file, output_file):
    """Reads the mpileup file and writes another file that is human-readable"""
    try:
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            for line in infile:
                generate_base = Base_Counter(line.strip())
                outfile.write(generate_base + "\n")
        print(f"{accession}: pileup created")
    except FileNotFoundError as ex:
        print(f"{ex}: File not found")
