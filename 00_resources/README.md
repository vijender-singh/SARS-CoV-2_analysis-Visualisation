

#### PLEASE NOTE:
1. *Scripts composed here are ment to be executed on Xanadu cluster (Bioinformatics HPC facility at University of Connecticut, running SLURM for resource management an job scheduling).  Hence the resource header is written for a SLURM scheduler.  Users can change this section to make it compatible to their own job scheduling system.*
2. *On Xanadu cluster the software/package/application are installed as modules which may be a different case on other systems*

### Reference Sequence and Annotation
The resource files used here were downloaded from NCBI's `SARS-COV-2` resource [webpage](https://www.ncbi.nlm.nih.gov/sars-cov-2/). The reference sequence `Wuhan-Hu-1_NC_045512.2.fasta` was obtained from [https://www.ncbi.nlm.nih.gov/nuccore/NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512).  The genome annotation `gff3`file is downloaded from [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz).

### Indexing
The reference sequence was `indexed` with mappers  `minimpa2` and `bwa-mem2`. A fasta index was also generated using `samtools faidx`. The script `refrence_index.sh` can be run to download annotation file and to create these indexes. 

### Obtaining the genome sequence
Go to the GenBank sequence (page)[https://www.ncbi.nlm.nih.gov/nuccore/NC_045512]. The displayed sequence is in GenBAnk format.  To download in `fasta` format select the `Send to` option displayed on the top right of the page below SEARCH bar. Select `Complete Record` and then `File` under Choose Destination and `FASTA` under Format and then click on `Create File`.  This will download the genome refernce file in fasta format.  It will download the file to your local computer.