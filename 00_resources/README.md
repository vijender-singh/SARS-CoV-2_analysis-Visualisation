

#### PLEASE NOTE:
*Scripts composed here are ment to be executed on Xanadu cluster (Bioinformatics HPC facility at University of Connecticut, running SLURM for resource management an job scheduling).  Hence the resource header written for a SLURM scheduler.  Users can change this section to make it compatible to their own job scheduling system*


### Reference Sequence and Annotation
The resource files used here are downloaded from NCBI's `SARS-COV-2` resource [webpage](https://www.ncbi.nlm.nih.gov/sars-cov-2/). The reference sequence `Wuhan-Hu-1_NC_045512.2.fasta` was obtained from [https://www.ncbi.nlm.nih.gov/nuccore/NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512).  The genome annotation `gff3`file is downloaded from [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz).

### Indexing
The reference sequence was `indexed` with mappers  `minimpa2` and `bwa-mem2`. A fasta index was also generated using `samtools faidx`. The script `refrence_index.sh` can be run to create these indexes. 