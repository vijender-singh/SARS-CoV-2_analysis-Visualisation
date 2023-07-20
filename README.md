## SARS-CoV-2_analysis-Visualisation

## Workflow Execution
> Copy the `config_file.txt` and make appropriate changes to it according to instructions described bellow in INPUT section. For each seperate analysis a new copy of the file can be made and be named as per users need.  Also the resource directory can be updated based on the instruction given in the `README.md` file in `00_resources` directory.  The file path to the resource directory can then be added to the config file.  Once the config file and resource directory are updated as required. Run the analysis as
> 
**`sh run.sh path/to/config/file/config_file.txt`**

## INPUT

The workflow runs in one of the 3 modes `PE`,`SE` and `PEtaSE` specified by the user.

> **PE**: In **`PE`** mode the input data is paired end sequecing with **`R1`** and **`R2`** files in **`.fastq.gz`** format. 

> **SE**: In **`SE`** mode the data will be treated as single-end sequencing even if the input data is paired-end.  All the files will be merged to generate a single fastq file.  

> **PEtaSE**: (**P**aired **E**nd **t**reated **a**s **S**ingle **E**nd ) This format would be ideal if the **`R1`** and **`R2`** reads of **`PE`** data overlap with each other.  In this format the reads will be merged (with minimum 6bps of overlap) to create a single contigous read.  This read will be treated as a **`SE`** data. The reads which failed to have an overlap (due to trimming) will be treated as single end reads.

In order to run the workflow please make adequate changes in the **`config_file.txt`** file.


1. **`DataDir="/path/to/data/directory"`**

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

* : Represent any strings present in fastq file names e.g. S1_L001 etc.  Keep the files compressed.

```

2. #### **`SampleNames=(Sample1 Sample2 Sample3 ... SampleN)`**
> Make sure the **`SampleNames`** matches the name of the directories containing the sample data as showed under `DataDir`.


3. #### **`SampleSources=(WW SI SI ... WW)`**
>
> Please indicate a Sample source for each sample.  The options are.
>
> **`WW`** : Waste Water Treatment.
> 
> **`SI`** : Sample Isolate/ clinical isolate.


4. #### **`AnalysisType=PEtaSE`**

> This parameter will indicate what kind of analysis has to be done.  The options are **`PE`**, **`SE`** and **`PEtaSE`**,  Each option is described above. 


5. #### **`MultiTypeAnalysis=`**
>
> Use `MultiTypeAnalysis` option only if each sample has to be processed in a different mode otherwise leave it blank. The syntax to specify it is 
> **`MultiTypeAnalysis=(PE SE PE ... PEtaSE)`**.


6. #### **`OutDir="/Path/to/output/Dir"`**
> **`OutDir`** directory will be created if that doesnot exist.


7. #### **`ProjectID="SKA12Jun2021"`**

> Give an ID to the project. 

8. #### **`ResourceDir=/path/to/resourrce/directory/00_resources`**

> Provide path to the resource directory (00_resources) included here.



## OUTPUT

The output will generate the following directory structure. The output is an example of running data with `example_config_file.txt`.  The data is analysed in `PEtaSE` mode.

```
Kendra_ISG_TEST             
├── Kendra_ISG_PEtaSE       
│   ├── MapDir-30Mar2021
│   ├── rawfastqcDir-30Mar2021
│   ├── trimDir-30Mar2021
│   ├── trimfastqcDir-30Mar2021
│   └── VarDir-30Mar2021
└── log-Kendra_ISG
    └── log-30_Mar_13hr-58min-37sec_PEtaSE
 ```
 
 ### Lets Look closely at the directory structure.
 
 ```
 Kendra_ISG_TEST
 
 This was the OutDir of the config_file.txt.  This Directory will contain 2 directories names based on the ProjectID (here: Kendra_ISG).
 (1): ProjectID_AnalysisType (Here : Kendra_ISG_PEtaSE).  This directory contain output from each individual steps in a seperate directory.
        (a) rawfastqcDir-dateofAnalysis : FastQC report of the raw files
        (b) trimDir-dateofAnalysis      : Trimmomatic trimming of the reads
        (c) trimfastqcDir-dateofAnalysis: FastQc of trimmed reads
        (d) VarDir-dateofAnalysis       : This directory contains, read depth information, bedgraph file of coverage, VCf files with variants and annotation of VCF file using the GFF file.  The annotation will have the information on gene the SNP is located in and also amino acid changes due to the SNP and more.
 (2): log-TimeStampOfAnalysis_AnalysisType:  ( here: log-30_Mar_13hr-58min-37sec_PEtaSE).  This contain all the log info of each run.  Each run will create a seperate log folder based on time stamp.
        
 ### Data Visualisation
Use the remote desktop app and connect to R-shiny app server on the VM 'sarcov2.cam.uchc.edu' (use CAM credentials to login) and add the `VarDir-dateofAnalysis` folder of the current analysis, and this will update the plots and graphs with new data included.
 ```
 
 
 

