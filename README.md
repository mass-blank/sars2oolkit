# sars2oolkit
Python toolkit to help identify SARS-CoV2 genome, map sequencer reads onto reference genome and more.


There are several programs that need to be installed for this script to work.

I recommend installing anaconda and creating a new python environment

conda config --add channels defaults

conda config --add channels bioconda

conda config --add channels conda-forge

conda config --set channel_priority strict

conda create -n bio

conda install biopython

conda install bowtie2 

conda install ... 


BINARY DEPENDENCIES

bowtie2

sra-toolkit

samtools

bcftools

fastqc 

fastv 

parallel-fastq-dump

entrez-direct


PYTHON DEPENDENCIES

numpy

pandas

biopython
