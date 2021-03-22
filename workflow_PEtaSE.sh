#!/bin/bash

source config_file.txt

SampleName=${1}
SampleSource=${2}
ProjectDir=${3}
OutDir=${4}
ProjectID=${5}

SampleDir=${ProjectDir}/${SampleName}t


queue=general
QOS=general
date

logTimeStamp=`date +%d"_"%b"_"%H"hr-"%M"min-"%S"sec"`
datestamp=`date +%d%b%Y`

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
        #SBATCH -J PEtaSE
        hostname
        date
        cd ${SampleDir}

        module load ea-utils/1.04.807

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
        	echo -e "Reads will be joined based on minimum overlap of 6bps"
        	#echo -e "${SampleName}_R1.fastq"
        	#echo -e "${SampleName}_R2.fastq"
        	gunzip *.fastq.gz

        	fastq-join *R1*.fastq *R2*.fastq -o ${SampleName}.%.fastq

        	cat ${SampleName}.*.fastq >> ${SampleName}.fastq

        else
        	echo -e "There is ${R1files} input file for R1 and ${R2files} for R2"
        	echo -e "R1 and R2 group of files will be merged and will be renamed as"
        	echo -e "${SampleName}_R1.fastq"
        	echo -e "${SampleName}_R2.fastq"
        	echo -e "Reads from the R1 and R2 files will be stitched based on minimum overlap of 6bps"
        	zcat `ls *R1*.fastq.gz | sort -` >> ${SampleName}_R1.fastq
        	zcat `ls *R2*.fastq.gz | sort -` >> ${SampleName}_R2.fastq

        	fastq-join ${SampleName}_R1.fastq ${SampleName}_R2.fastq -o ${SampleName}.%.fastq

        	cat ${SampleName}.*.fastq >> ${SampleName}_PEtaSE.fastq
        fi

        	readsinital=`wc -l ${SampleName}_R1.fastq`
        	readsmerged=`wc -l ${SampleName}.join.fastq`
        	prcnt=$(( 100*${readsmerged}/${readsinital} ))

        	echo -e "       \t Total Reads \t reads_joined "
        	echo -e " ========================================================"
        	echo -e "counts \t ${readsinital} \t ${readsmerged}"
        	echo -e "percent\t                \t ${prcnt}%"


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

        fastqc -o ${RawfastqcDir} --threads 4  ${SampleName}_PEtaSE.fastq

		module rm fastqc/0.11.7
		module load Trimmomatic/0.39

		java -jar $Trimmomatic PE -threads 4 \
		    ${SampleName}_PEtaSE.fastq \
		    ${trimDir}/trim_${SampleName}_PEtaSE.fastq\
		    ILLUMINACLIP:${resourceDir}/TruSeq3-SE.fa:2:30:10 \
		    SLIDINGWINDOW:4:25 \
		    MINLEN:45

		module rm Trimmomatic/0.39
        module load fastqc/0.11.7

        fastqc -o ${TrimfastqcDir} --threads 4  ${trimDir}/trim_${SampleName}_PEtaSE.fastq

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
				sr ${resourceDir}/minimap2/Wuhan-Hu-1_NC_045512.2.mmi  ${trimDir}/trim_${SampleName}_PEtaSE.fastq \
				| samtools view -h -F 4 | samtools sort --threads 4 -o ${MapDir}/minimap_${SampleName}_PEtaSE.bam

			cd ${MapDir}

			samtools index minimap_${SampleName}_PEtaSE.bam

			samtools depth -d 0 -Q 30 -q 0 -aa minimap_${SampleName}_PEtaSE.bam > ${OutDir}/Depth_minimap_${SampleName}_PEtaSE.txt

			bedtools genomecov -bg -ibam minimap_${SampleName}_PEtaSE.bam >> ${OutDir}/${SampleName}_PEtaSE.bdg

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
	        	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -C 5 --pooled-continuous ${MapDir}/minimap_${SampleName}_PEtaSE.bam > ${SampleName}_PEtaSE.vcf
	        else
	        	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -F 0.2 ${MapDir}/minimap_${SampleName}_PEtaSE.bam > ${SampleName}_PEtaSE.vcf
	        fi

	        #snpEff : Wuhan-Hu-1_NC_045512.2

	        module load snpEff/4.3q

	        java -Xmx20g -jar $SNPEFF Wuhan-Hu-1_NC_045512.2 ${SampleName}_PEtaSE.vcf > ${SampleName}__PEtaSE_Anno.vcf


	        echo "$0 $@"
	
	VARCALL`






















