
OutDir="Path/to/outputDIr"
ProjectDir="Path/to/ProjectDir"
SampleName="SampleName"

SampleSource="SampleSourceType WW/SI"
ReadType="PE/SE/PEtaSE"
ProjectID="PID"
ReRunID="18Mar2021"

resourceDir="path/to/00_resources"



queue=general
QOS=general
date

logTimeStamp=`date +%d"_"%b"_"%H"hr-"%M"min-"%S"sec"`
datestamp=`date +%d%b%Y`

SampleDir=${ProjectDir}/data/${SampleName}
ProjectLogDir=${OutDir}/log-${ProjectID}
RawfastqcDir=${OutDir}/rawfastqcDir-${datestamp}
TrimfastqcDir=${OutDir}/trimfastqcDir-${datestamp}
trimDir=${OutDir}/trimDir-${datestamp}
MapDir=${OutDir}/MapDir-${datestamp}
VarDir=${OutDir}VarDir-${datestamp}

if [ -z "${ReRunID}" ];then
	mkdir -p ${rawfastqcDir} ${trimDir} ${MapDir} ${VarDir}
	chmod 777 ${rawfastqcDir} ${trimDir} ${MapDir} ${VarDir}
fi

LogDir=${ProjectLogDir}/log-${logTimeStamp}

if [ ! -d "${LogDir}" ]; then
    mkdir -p "${LogDir}"
    chmod 777 "${LogDir}"
fi

### Merging Files
cd ${SampleDir}

jid1=`sbatch <<- MERGE | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p ${queue}
        #SBATCH -q ${QOS}
        #SBATCH -c 1
        #SBATCH --mem=5G
        #SBATCH -o ${LogDir}/merge-%j.out
        #SBATCH -e $debugdir/merge-%j.err
        #SBATCH -J PE_merge_files
        hostname
        date

        cd ${SampleDir}

        R1files=`ls *R1*.fastq.gz 2>/dev/null | wc -l`
        R2files=`ls *R2*.fastq.gz 2>/dev/null | wc -l`

        if [ ${R1files} != ${R2files} ]
        then
        	echo -e "The number of input R1 and R2 files for the sample are not same. "
        	echo -e "Please check the input files"
        	echo -e "Aborting at MERGE step"
        	exit 1
        elif [ ${R1files} -eq 1 ] && [ ${R2files} -eq 1 ]
        then
        	echo -e "There is 1 input file each for R1 and R2"
        	echo -e "Files will not be merged but will be renamed as"
        	echo -e "${SampleName}_R1.fastq"
        	echo -e "${SampleName}_R2.fastq"
        	zcat `ls *R1*.fastq.gz ` >> ${SampleName}_R1.fastq
        	zcat `ls *R2*.fastq.gz ` >> ${SampleName}_R2.fastq
        else
        	echo -e "There is ${R1files} input file for R1 and ${R2files} for R2"
        	echo -e "R1 and R2 group of files will be merged and will be renamed as"
        	echo -e "${SampleName}_R1.fastq"
        	echo -e "${SampleName}_R2.fastq"
        	zcat `ls *R1*.fastq.gz ` >> ${SampleName}_R1.fastq
        	zcat `ls *R2*.fastq.gz ` >> ${SampleName}_R2.fastq
    	    zcat `ls *R1*.fastq.gz | sort -` >> ${SampleName}_R1.fastq
        	zcat `ls *R2*.fastq.gz | sort -` >> ${SampleName}_R2.fastq
        fi

       	reads1=`wc -l < ${SampleName}_R1.fastq`
		reads2=`wc -l < ${SampleName}_R2.fastq`
		total_reads=$(( $reads1 / 4 ))
		if [ $reads1 == $reads2 ];then
		          echo "File merging Authenticated"
		          echo "Total reads are ${total_reads} in R1 and R2 files"
		else
		          echo "There was issue merging the files, The number of reads in both files donot match...!!!!"
		          exit 1
		fi



		        echo "$0 $@"
		MERGE`


depend1="afterok:${jid1}"

jid2=`sbatch <<- TRIMMOTRIM | egrep -o -e "\b[0-9]+$"
	        #!/bin/bash -l
	        #SBATCH -p ${queue}
	        #SBATCH -q ${QOS}
	        #SBATCH -c 1
	        #SBATCH --mem=20G
	        #SBATCH -d ${depend1}
	        #SBATCH -o ${LogDir}/QC_TRIM_QC-%j.out
	        #SBATCH -e $debugdir/QC_TRIM_QC-%j.err
	        #SBATCH -J QC_TRIM_QC
	        hostname
	        date

	        cd ${SampleDir}

	        module load fastqc/0.11.7

	        fastqc -o ${RawfastqcDir} --threads 4  ${SampleName}_R1.fastq ${SampleName}_R2.fastq

			module rm fastqc/0.11.7
			module load Trimmomatic/0.39

			java -jar $Trimmomatic PE -threads 4 \
			    ./${SampleName}_R1.fastq \
			    ./${SampleName}_R2.fastq \
			    ${trimDir}/trim_${SAM}_R1.fastq.gz ${trimDir}/singles_trim_${SAM}_R1.fastq \
			    ${trimDir}/trim_${SAM}_R2.fastq.gz ${trimDir}/singles_trim_${SAM}_R2.fastq \
			    ILLUMINACLIP:${resourceDir}/TruSeq3-PE-2.fa:2:30:10 \
			    SLIDINGWINDOW:4:25 \
			    MINLEN:45

			module rm Trimmomatic/0.39
	        module load fastqc/0.11.7

	        fastqc -o ${TrimfastqcDir} --threads 4  ${trimDir}/trim_${SampleName}_R1.fastq ${trimDir}/trim_${SampleName}_R2.fastq

			module rm fastqc/0.11.7
			
			gzip *.fastq

	        echo "$0 $@"
	TRIMMOTRIM`


depend2="afterok:${jid2}:${jid1}"

jid3=`sbatch <<- MAPPING | egrep -o -e "\b[0-9]+$"
	        #!/bin/bash -l
	        #SBATCH -p ${queue}
	        #SBATCH -q ${QOS}
	        #SBATCH -c 4
	        #SBATCH --mem=20G
	        #SBATCH -d ${depend2}
	        #SBATCH -o ${LogDir}/map_${SampleName}-%j.out
	        #SBATCH -e ${LogDir}/map_${SampleName}-%j.err
	        #SBATCH -J map_${SampleName}
	        hostname
	        date

	        module load minimap2/2.17
			module load samtools
			module load bedtools/2.29.0

			minimap2 -t 4 -a -x \
				sr ${resourceDir}/minimap2/Wuhan-Hu-1_NC_045512.2.mmi  ${trimDir}/trim_${SampleName}_R1.fastq ${trimDir}/trim_${SampleName}_R2.fastq \
				| samtools view -h -F 4 | samtools sort --threads 4 -o ${MapDir}/minimap_${SampleName}_PE.bam

			cd ${MapDir}

			samtools index minimap_${SampleName}_PE.bam

			samtools depth -d 0 -Q 30 -q 0 -aa minimap_${SampleName}_PE.bam > ${OutDir}/Depth_minimap_${SampleName}_PE.txt

			bedtools genomecov -bg -ibam minimap_${SampleName}_PE.bam >> ${OutDir}/${SampleName}_PE.bdg

	        echo "$0 $@"
	
	MAPPING`

depend3="afterok:${jid3}"


jid4=`sbatch <<- VARCALL | egrep -o -e "\b[0-9]+$"
	        #!/bin/bash -l
	        #SBATCH -p ${queue}
	        #SBATCH -q ${QOS}
	        #SBATCH -c 4
	        #SBATCH --mem=25G
	        #SBATCH -d ${depend3}
	        #SBATCH -o ${LogDir}/VarCall_${SampleName}-%j.out
	        #SBATCH -e ${LogDir}/VarCall_${SampleName}-%j.err
	        #SBATCH -J VarCall_${SampleName}
	        hostname
	        date

	        cd ${VarDir}

	        sampleSource=${sampleSource}


	        module load freebayes/1.3.4
	        if [ ${sampleSource} = "WW" ]; then
	        	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -C 5 --pooled-continuous ${MapDir}/minimap_${SampleName}_PE.bam > ${SampleName}_PE.vcf
	        else
	        	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -F 0.2 ${MapDir}/minimap_${SampleName}_PE.bam > ${SampleName}_PE.vcf
	        fi

	        #snpEff : Wuhan-Hu-1_NC_045512.2

	        module load snpEff/4.3q

	        java -Xmx20g -jar $SNPEFF Wuhan-Hu-1_NC_045512.2 ${SampleName}_PE.vcf > ${SampleName}_PE_Anno.vcf


	        echo "$0 $@"
	
	VARCALL`



















