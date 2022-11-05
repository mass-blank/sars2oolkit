from pathlib import Path


class Accession:
    """Creates an Accession object with paths to file names."""
    def __init__(self, accession: str):
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


class Line:
    def __init__(self, nucleotide, bigA, bigC, bigG, bigT):
        self.nucleotide = int(nucleotide)
        self.A: int = int(bigA)
        self.C: int = int(bigC)
        self.G: int = int(bigG)
        self.T: int = int(bigT)
