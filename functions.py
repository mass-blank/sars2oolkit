import os
import shutil
import subprocess
import json


def call_mutations(accession):
    args = f"bcftools mpileup -d 150000 -Ou -f SARS-CoV-2_reference.fasta {accession}.sorted.bam | bcftools call -mv -Ob -o {accession}.bcf"
    subprocess.run(args, shell=True)
    # args2 = "bcftools view - i '%QUAL>=15' - H calls.bcf
    # process all alleles no file: bcftools mpileup  -a AD -Ou -f SARS-CoV-2_reference.fasta SRR12019375.sorted.bam | bcftools call -mA


def gen_pileup(accession, nt_start=None, nt_stop=None):
    args = f"samtools mpileup -d 150000 -r NC_045512.2:{nt_start}-{nt_stop} -o {accession}_pileup.txt {accession}.sorted.bam"
    subprocess.run(args, shell=True)


def view_mutations(accession):
    args = "bcftools query -f '" + accession + r" %REF%POS%ALT\n' " + accession + ".bcf"
    return subprocess.run(args, shell=True, capture_output=True, text=True)


def isPairedSRA(accession):
    filename = os.path.abspath(accession)
    try:
        contents = subprocess.check_output(["fastq-dump", "-X", "1", "-Z", "--split-spot", filename])
        if (contents.count(b'\n') == 4):
            return False
        elif (contents.count(b'\n') == 8):
            return True
        else:
            pass
    except subprocess.CalledProcessError:
        raise Exception("Error running fastq-dump on ", accession)


def bow_tie(accession, sra_number_1=None, sra_number_2=None):
    if isPairedSRA(accession):
        args = "bowtie2 " + " -x " + " bowtie " + " -1 " + str(sra_number_1) + " -2 " + str(sra_number_2) + " -S " + accession + ".sam"
        subprocess.run(args, shell=True)
    else:
        args = "bowtie2 " + " -x " + " bowtie " + " -U" + accession + ".fastq.gz" + " -S " + accession + ".sam"
        subprocess.run(args, shell=True)


def sam_tools_view(accession):
    args = "samtools " + "view " + "-S " + "-b " + accession + ".sam " + "> " + accession + ".bam"
    subprocess.run(args, shell=True)


def sam_tools_sort(accession):
    args = "samtools " + "sort " + accession + ".bam " + "-o " + accession + ".sorted.bam"
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
    subprocess.run(["prefetch", accession])


def fastqc_func(accession, fastqfile_1=None, fastqfile_2=None):
    if isPairedSRA(accession):
        args = f"fastqc {fastqfile_1} {fastqfile_2}"
        subprocess.run(args, shell=True)
    else:
        args = f"fastqc {accession}.fastq.gz"
        subprocess.run(args, shell=True)


def fastv_func(accession, fastqfile_1=None, fastqfile_2=None):
    if isPairedSRA(accession):
        args = "fastv" + " --in1 " + str(fastqfile_1) + " --in2 " + str(fastqfile_2) + " -g " + "SARS-CoV-2.genomes.fa" + " -k " + "SARS-CoV-2.kmer.fa" + " -h " + accession + ".html" + " -j " + accession + ".json"
        subprocess.run(args, shell=True)
    else:
        args = "fastv" + " -i " + accession + ".fastq.gz" + " -g " + "SARS-CoV-2.genomes.fa" + " -k " + "SARS-CoV-2.kmer.fa" + " -h " + accession + ".html" + " -j " + accession + ".json"
        subprocess.run(args, shell=True)


def check_positive(accession):
    with open(accession + ".json", "r") as file:
        data = json.loads(file.read())
        return(data["kmer_detection_result"]["result"])


def split_range(nt_range):
    return nt_range.split('-')


def is_full():
    total, used, free = shutil.disk_usage('/')
    free = free / 2**20
    if free < 900:
        print("There is only %f free on the hard disk. Deleting..." % free)
        return True
    else:
        return False
