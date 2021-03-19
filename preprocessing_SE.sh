
OutDir="Path/to/outputDIr"
ProjectDir="Path/to/ProjectDir"
SampleName="SampleName"
SampleDir=${ProjectDir}/${SampleName}
SampleSource="SampleSourceType WW/SI"
ReadType="PE/SE/PEtaSE"
ProjectID="PID"
ReRunID="18Mar2021"



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
        #SBATCH -J SE_merge
        hostname
        date
        cd ${SampleDir}

        Readfiles=`ls *.fastq.gz 2>/dev/null | wc -l`
        	echo -e "There is/are ${Readfiles} input fastq.gz files"
        	echo -e "Files will be merged to a single file"
        	zcat `ls *.fastq.gz ` >> ${SampleName}_SE.fastq


        echo "$0 $@"
MERGE`



depend1="afterok:${jid1}"


jid2=`sbatch <<- TRIMMOTRIM | egrep -o -e "\b[0-9]+$"
        #!/bin/bash -l
        #SBATCH -p ${queue}
        #SBATCH -q ${QOS}
        #SBATCH -c 4#SBATCH --mem=20G
        #SBATCH -d ${depend1}
        #SBATCH -o ${LogDir}/QC_TRIM_QC-%j.out
        #SBATCH -e $debugdir/QC_TRIM_QC-%j.err
        #SBATCH -J QC_TRIM_QC
        hostname
        date
        cd ${SampleDir}

	    module load fastqc/0.11.7

        fastqc -o ${RawfastqcDir} --threads 4  ${SampleName}_SE.fastq

		module rm fastqc/0.11.7
		module load Trimmomatic/0.39

		java -jar $Trimmomatic PE -threads 4 \
		    ${SampleName}_SE.fastq \
		    ${trimDir}/trim_${SampleName}_SE.fastq\
		    ILLUMINACLIP:${resourceDir}/TruSeq3-SE.fa:2:30:10 \
		    SLIDINGWINDOW:4:25 \
		    MINLEN:45

		module rm Trimmomatic/0.39
        module load fastqc/0.11.7

        fastqc -o ${TrimfastqcDir} --threads 4  ${trimDir}/trim_${SampleName}_SE.fastq

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
				sr ${resourceDir}/minimap2/Wuhan-Hu-1_NC_045512.2.mmi  ${trimDir}/trim_${SampleName}_SE.fastq \
				| samtools view -h -F 4 | samtools sort --threads 4 -o ${MapDir}/minimap_${SampleName}_SE.bam

			cd ${MapDir}

			samtools index minimap_${SampleName}_SE.bam

			samtools depth -d 0 -Q 30 -q 0 -aa minimap_${SampleName}_SE.bam > ${OutDir}/Depth_minimap_${SampleName}_SE.txt

			bedtools genomecov -bg -ibam minimap_${SampleName}_SE.bam >> ${OutDir}/${SampleName}_SE.bdg

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
	        	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -C 5 --pooled-continuous ${MapDir}/minimap_${SampleName}_SE.bam > ${SampleName}_SE.vcf
	        else
	        	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -F 0.2 ${MapDir}/minimap_${SampleName}_SE.bam > ${SampleName}_SE.vcf
	        fi

	        #snpEff : Wuhan-Hu-1_NC_045512.2

	        module load snpEff/4.3q

	        java -Xmx20g -jar $SNPEFF Wuhan-Hu-1_NC_045512.2 ${SampleName}_SE.vcf > ${SampleName}_SE_Anno.vcf


	        echo "$0 $@"
	
	VARCALL`























