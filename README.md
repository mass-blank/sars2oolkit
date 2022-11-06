# sars2oolkit
###### Python toolkit to help identify SARS-CoV2 genome, map sequencer reads onto reference genome and more.


There are several programs that need to be installed for this script to work.

I recommend installing anaconda and creating a new python environment
````
conda config --add channels defaults

conda config --add channels bioconda

conda config --add channels conda-forge

conda config --set channel_priority strict

conda create -n bio

conda install biopython

conda install bowtie2 

conda install ... 
````


### BINARY DEPENDENCIES

1. [x] bowtie2
2. [x] bwa
3. [x] sra-toolkit
4. [x] samtools
5. [x] bcftools
6. [x] fastqc
7. [x] fastv
8. [x] parallel-fastq-dump
9. [x] entrez-direct

### PYTHON DEPENDENCIES
1. [x] numpy
2. [x] pandas
3. [x] biopython

### How to run the program?

#### Download Accessions
```commandline
python checksars2.py -i [accession_file.txt] -d
```
Run the above command several times until all files have been processed and downloaded.
#### Align Accessions 
```commandline
python checksars2.py -i [accession_file.txt] -b
```