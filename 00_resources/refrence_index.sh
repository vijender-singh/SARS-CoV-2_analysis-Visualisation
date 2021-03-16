#!/bin/bash
#SBATCH -J Indexing
#SBATCH -p general
#SBATCH -q general
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -o ./index-%j.out


hostname # Compute-node information is collecated on which the script will be executed.
date

mkdir bwamem2 minimap2 annotation

###############################
# Annotation Download
###############################

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz

mv GCF_009858895.2_ASM985889v3_genomic.gff.gz ./annotation/Wuhan-Hu-1_NC_045512.2.gff.gz 


RefFile=./Wuhan-Hu-1_NC_045512.2.fasta

#create Output Directories

###############################
# minimap2 indexing
###############################

module load minimap2/2.17
minimap2 -d ./minimap2/Wuhan-Hu-1_NC_045512.2.mmi ${RefFile}
module rm minimap2/2.17

###############################
# bwa-mem2 indexing
###############################

module load bwa-mem2/2.1
bwa-mem2 index -p ./bwamem2/Wuhan-Hu-1_NC_045512.2 ${RefFile}
module rm bwa-mem2/2.1

###############################
# samtool indexing
###############################

module load samtools/1.10
samtools faidx ${RefFile}


date