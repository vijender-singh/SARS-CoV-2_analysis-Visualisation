## SARS-CoV-2_analysis-Visualisation

The workflow runs in 3 modes `PE`,`SE` and `PEtaSE`.

**PE**: In **`PE`** mode the input data is paired end sequecing with **`R1`** and **`R2`** files in **`.fastq.gz`** format. 

**SE**: In **`SE`** mode the data will be treated as single-end sequencing even if the input data is paired-end.  All the files will be merged to generate a single fastq file.  

**PEtaSE**: (**P**aired **E**nd **t**reated **a**s **S**ingle **E**nd ) This format would be ideal if the **`R1`** and **`R2`** reads of **`PE`** data overlap with each other.  In this format the reads will be merged (with minimum 6bps of overlap) to create a single contigous read.  This read will be treated as a **`SE`** data. The reads which failed to have an overlap (due to trimming) will be treated as single end reads.

In order to run the workflow please maked adeaquate changes in the **`config_file.txt`** file.
