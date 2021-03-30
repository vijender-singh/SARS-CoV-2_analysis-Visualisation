## SARS-CoV-2_analysis-Visualisation

The workflow runs in 3 modes `PE`,`SE` and `PEtaSE`.

**PE**: In **`PE`** mode the input data is paired end sequecing with **`R1`** and **`R2`** files in **`.fastq.gz`** format. 

**SE**: In **`SE`** mode the data will be treated as single-end sequencing even if the input data is paired-end.  All the files will be merged to generate a single fastq file.  

**PEtaSE**: (**P**aired **E**nd **t**reated **a**s **S**ingle **E**nd ) This format would be ideal if the **`R1`** and **`R2`** reads of **`PE`** data overlap with each other.  In this format the reads will be merged (with minimum 6bps of overlap) to create a single contigous read.  This read will be treated as a **`SE`** data. The reads which failed to have an overlap (due to trimming) will be treated as single end reads.

In order to run the workflow please maked adequate changes in the **`config_file.txt`** file.

**`DataDir="/path/to/data/directory"`**

*The data directory should have following heirarchy *
```
DataDir/
├── Sample1
│   ├── Sample1*_R1.fastq.gz
│   └── Sample1*_R2.fastq.gz
├── Sample2
│   ├── Sample2*_R1.fastq.gz
│   └── Sample2*_R2.fastq.gz
├── Sample3
│   ├── Sample3*_R1.fastq.gz
│   └── Sample3*_R2.fastq.gz
.
.
.
└── SampleN
    ├── SampleN*_R1.fastq.gz
    └── SampleN*_R2.fastq.gz

* : Represent any strings present is sample names e.g. S1_L001 etc.
```

