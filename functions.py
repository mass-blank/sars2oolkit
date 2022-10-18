import json
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


def xml_parse(accession, search):
    args = f"esummary -db sra -id {accession}"
    myxml = subprocess.run(args, shell=True, capture_output=True, text=True)
    myroot = ET.fromstring(myxml.stdout)
    result = [result.text for result in myroot.findall(f'.//{search}')]
    return result


def call_mutations(accession):
    args = f"bcftools mpileup -d 150000 -Ou -f SARS-CoV-2_reference.fasta {accession}.sorted.bam | bcftools call -mv -Ob -o {accession}.bcf"
    subprocess.run(args, shell=True)


def gen_pileup(accession):
    args = (
        f"samtools mpileup -d 150000 -o {accession}_pileup.txt {accession}.sorted.bam"
    )
    try:
        subprocess.run(args, shell=True)
    except subprocess.CalledProcessError as ex:
        print(f"{ex.output} Couldn't run mpileup")


def view_mutations(accession):
    args = "bcftools query -f '" + accession + \
        r" %REF%POS%ALT\n' " + accession + ".bcf"
    return subprocess.run(args, shell=True, capture_output=True, text=True)


def fastq_exists(accession):
    try:
        result = subprocess.run(
            ["prefetch", "--type", "fastq", accession], capture_output=True
        )
        re_result = re.findall(r"'.*[gz]'", result.stderr.decode("utf-8"))
        re_result = [i.strip("'") for i in re_result]
        re_result = list(dict.fromkeys(re_result))
        current_dir = os.path.abspath(os.getcwd())
        if len(re_result) == 1:
            infile = Path(accession + "/" + re_result[0])
            outfile = Path(accession + ".fastq.gz")
            os.rename(infile, os.path.join(current_dir, outfile))
            print(f"{outfile.name} has finished downloading")
        elif len(re_result) == 2:
            infile1 = Path(accession + "/" + re_result[0])
            infile2 = Path(accession + "/" + re_result[1])
            outfile1 = Path(accession + "_1.fastq.gz")
            outfile2 = Path(accession + "_2.fastq.gz")
            os.rename(infile1, os.path.join(current_dir, outfile1))
            os.rename(infile2, os.path.join(current_dir, outfile2))
            print(f"{outfile1.name} and {outfile2.name} have finished downloading")
        elif not re_result:
            subprocess.run(["prefetch", accession])
        else:
            print("Fastq_exists else conditional")
    except subprocess.CalledProcessError:
        raise Exception("Error running fastq-dump on ", accession)


def isPairedSRA(accession):
    filename = Path(os.path.abspath(accession))
    if filename.is_file():
        contents = subprocess.check_output(
            ["fastq-dump", "-X", "1", "-Z", "--split-spot", filename]
        )
        if contents.count(b"\n") == 4:
            return False
        elif contents.count(b"\n") == 8:
            return True
        else:
            pass
    else:
        print(filename.name + " does not exist.")


def bow_tie(accession, fastq_file_1=None, fastq_file_2=None, fastq_file=None):
    if fastq_file_1 is not None and fastq_file_2 is not None:
        if fastq_file_1.is_file() and fastq_file_2.is_file():
            args = f"bowtie2 -p 16 -x bowtie -1 {str(fastq_file_1)} -2 {str(fastq_file_2)} -S {accession}.sam"
            subprocess.run(args, shell=True)
        elif fastq_file.is_file():
            args = f"bowtie2 -p 16 -x bowtie -U {accession}.fastq.gz -S {accession}.sam"
            subprocess.run(args, shell=True)
        else:
            print("No FASTQ files.")
    else:
        print("bow_tie arguments are none")


def sam_tools_view(accession):
    args = f"samtools view -S -b {accession}.sam > {accession}.bam"
    subprocess.run(args, shell=True)


def sam_tools_sort(accession):
    args = f"samtools sort {accession}.bam -o {accession}.sorted.bam"
    subprocess.run(args, shell=True)


def sam_tools_index(accession):
    args = "samtools " + "index " + accession + ".sorted.bam"
    subprocess.run(args, shell=True)


def fastq_func(accession):
    if isPairedSRA(accession):
        subprocess.run(["fastq-dump", "--gzip", accession])
    else:
        subprocess.run(["fastq-dump", "--gzip", "--split-e", accession])


def fetch_func(accession):
    subprocess.run(["prefetch", accession], stdout=subprocess.PIPE)


def fastqc_func(accession, fastq_file_1=None, fastq_file_2=None):
    if isPairedSRA(accession):
        args = f"fastqc {fastq_file_1} {fastq_file_2}"
        subprocess.run(args, shell=True)
    else:
        args = f"fastqc {accession}.fastq.gz"
        subprocess.run(args, shell=True)


def fastv_func(accession, fastq_file_1=None, fastq_file_2=None):
    if fastq_file_1 is not None and fastq_file_2 is not None:
        if fastq_file_1.is_file() and fastq_file_2.is_file():
            args = (
                "fastv"
                + " --in1 "
                + str(fastq_file_1)
                + " --in2 "
                + str(fastq_file_2)
                + " -g "
                + "SARS-CoV-2.genomes.fa"
                + " -k "
                + "SARS-CoV-2.kmer.fa"
                + " -h "
                + accession
                + ".html"
                + " -j "
                + accession
                + ".json"
            )
            subprocess.run(args, shell=True)
    else:
        args = (
            "fastv"
            + " -i "
            + accession
            + ".fastq.gz"
            + " -g "
            + "SARS-CoV-2.genomes.fa"
            + " -k "
            + "SARS-CoV-2.kmer.fa"
            + " -h "
            + accession
            + ".html"
            + " -j "
            + accession
            + ".json"
        )
        subprocess.run(args, shell=True)


def check_positive(json_file):
    with open(json_file, "r") as file:
        data = json.loads(file.read())
        return data["kmer_detection_result"]["result"]


def delete_accession(file, accession):
    with open(file, "r+") as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            if accession not in line.strip("\n"):
                file.write(line)
        file.truncate()


def mean_depth(json_file):
    with open(json_file, "r") as file:
        data = json.loads(file.read())
        return float(data["kmer_detection_result"]["mean_coverage"])


def is_full():
    total, used, free = shutil.disk_usage("/")
    free = free / 2**20
    if free < 900:
        print("There is only %f free on the hard disk. Deleting..." % free)
        return True
    else:
        return False
