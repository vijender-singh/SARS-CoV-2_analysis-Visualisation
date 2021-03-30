#!/bin/bash
#SBATCH -p general
#SBATCH -q general
#SBATCH --mem=2G
#SBATCH -J MScript


SampleName=${1}
SampleSource=${2}
DataDir=${3}
OutDir=${4}
ProjectID=${5}

SampleDir=${DataDir}/${SampleName}

resourceDir=${6}
logTimeStamp=${7}
LogDir=${8}
ProjectLogDir=${9}


queue="general"
qos="general"


date

datestamp=`date +%d%b%Y`
SampleDir=${DataDir}/${SampleName}
RawfastqcDir=${OutDir}/rawfastqcDir-${datestamp}
TrimfastqcDir=${OutDir}/trimfastqcDir-${datestamp}
trimDir=${OutDir}/trimDir-${datestamp}
MapDir=${OutDir}/MapDir-${datestamp}
VarDir=${OutDir}/VarDir-${datestamp}

if [ -z "${ReRunID}" ];then
	mkdir -p ${RawfastqcDir} ${trimDir} ${MapDir} ${VarDir} ${TrimfastqcDir}
	chmod 777 ${RawfastqcDir} ${trimDir} ${MapDir} ${VarDir} ${TrimfastqcDir}
fi

Statuslogfile=${LogDir}/Statuslogfile-${logTimeStamp}

#if [ ! -d "${LogDir}" ]; then
#    mkdir -p "${LogDir}"
#    chmod 777 "${LogDir}"
#fi

if [ -f ${LogDir}/errorMERGE ]
then
    rm ${LogDir}/errorMERGE
fi




### Merging Files
cd ${SampleDir}
mkdir -p PE

fwd="$(ls -1 *R1*.fastq.gz | wc -l)"
rev="$(ls -1 *R2*.fastq.gz | wc -l)"

if [ -f ./PE/${SampleName}_R1.fastq ]; then
	rm ./PE/${SampleName}_R?.fastq
fi


jid1=`sbatch --parsable -W  <<- MERGEALL
#!/bin/bash -l 
#SBATCH -p ${queue}
#SBATCH -q ${qos}
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH -o ${LogDir}/${SampleName}-merge-PE_%j.out
#SBATCH -J PE_merge_files
hostname
date

cd ${SampleDir}
echo ${fwd} ${rev}
if [ ${fwd} -ne ${rev} ];
then
	echo -e "The number of input R1 ${fwd} and R2 ${rev} files for the sample are not same. "
	echo -e "Please check the input files"
	echo -e "Aborting at MERGE step"
	exit 1
	touch ${LogDir}/errorMERGE
elif [ ${fwd} -eq 1 ] && [ ${rev} -eq 1 ];
then
	echo -e "There is 1 input file each for R1 and R2"
	echo -e "Files will not be merged but will be renamed as"
	echo -e "${SampleName}_R1.fastq"
	echo -e "${SampleName}_R2.fastq"
	zcat $( ls *R1*.fastq.gz ) >> PE/${SampleName}_R1.fastq
	zcat $( ls *R2*.fastq.gz ) >> PE/${SampleName}_R2.fastq
else
	echo -e "There is ${fwd} input file for R1 and ${rev} for R2"
	echo -e "R1 and R2 group of files will be merged and will be renamed as"
	echo -e "${SampleName}_R1.fastq"
	echo -e "${SampleName}_R2.fastq"
	zcat $( ls *R1*.fastq.gz | sort - ) >> PE/${SampleName}_R1.fastq
	zcat $( ls *R2*.fastq.gz | sort - ) >> PE/${SampleName}_R2.fastq
fi

	echo "$0 $@"
touch ${Statuslogfile}

MERGEALL`

echo ${jid1}

depend1="afterok:${jid1}"

echo "Start"

i=0; j=0
while [ ${i} -lt 1 ]
do
if [ -f ${Statuslogfile} ];then
        ((++i))
else
        ((++j))
	sleep 5
	echo $j
fi
done

echo "Ended"

cd ${SampleDir}

reads1=$(wc -l < PE/${SampleName}_R1.fastq)
reads2=$(wc -l < PE/${SampleName}_R2.fastq)

total_reads=$(($reads1 / 4))
if [ ${reads1} -eq ${reads2} ];then
        echo -e "Reads STATS : ${reads1}  ${reads2}  ${total_reads}"
        echo -e "File merging Authenticated"
        echo -e "Total reads are ${total_reads} in R1 and R2 files"
else
        echo -e "There was issue merging the files, The number of reads in both files donot match...!!!!"
        exit 1
	touch ${Statuslogfile}_ERROR
fi



jid2=`sbatch --parsable <<- TRIMMOTRIM 
#!/bin/bash -l
#SBATCH -p ${queue}
#SBATCH -q ${qos}
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -d ${depend1}
#SBATCH -o ${LogDir}/${SampleName}-QC_TRIM_QC-PE_%j.out
#SBATCH -J QC_TRIM_QC
hostname
date
cd ${SampleDir}
module load fastqc/0.11.7

fastqc -o ${RawfastqcDir} --threads 4  PE/${SampleName}_R1.fastq PE/${SampleName}_R2.fastq

module rm fastqc/0.11.7
module load Trimmomatic/0.39

java -jar /isg/shared/apps/Trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 4 \
	./PE/${SampleName}_R1.fastq \
	./PE/${SampleName}_R2.fastq \
	${trimDir}/trim_${SampleName}_R1.fastq ${trimDir}/singles_trim_${SampleName}_R1.fastq \
	${trimDir}/trim_${SampleName}_R2.fastq ${trimDir}/singles_trim_${SampleName}_R2.fastq \
	ILLUMINACLIP:${resourceDir}/TruSeq3-PE-2.fa:2:30:10 \
	SLIDINGWINDOW:4:25 \
	MINLEN:45

module rm Trimmomatic/0.39
module load fastqc/0.11.7

fastqc -o ${TrimfastqcDir} --threads 4  ${trimDir}/trim_${SampleName}_R1.fastq ${trimDir}/trim_${SampleName}_R2.fastq

module rm fastqc/0.11.7
			
#gzip *.fastq

echo "$0 $@"

TRIMMOTRIM`


depend2="afterok:${jid2}"

echo ${depend2}

jid3=`sbatch --parsable <<- MAPPING
#!/bin/bash -l	        
#SBATCH -p ${queue}
#SBATCH -q ${qos}
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -d ${depend2}
#SBATCH -o ${LogDir}/${SampleName}-map-PE_%j.out
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

samtools depth -d 0 -Q 30 -q 0 -aa minimap_${SampleName}_PE.bam > ${VarDir}/Depth_minimap_${SampleName}_PE.txt

bedtools genomecov -bg -ibam minimap_${SampleName}_PE.bam >> ${VarDir}/${SampleName}_PE.bdg

echo "$0 $@"
	
MAPPING`

depend3="afterok:${jid3}"

echo ${depend3}

jid4=`sbatch --parsable <<- VARCALL
#!/bin/bash -l
#SBATCH -p ${queue}
#SBATCH -q ${qos}
#SBATCH -c 4
#SBATCH --mem=25G
#SBATCH -d ${depend3}
#SBATCH -o ${LogDir}/${SampleName}-VarCall-PE_%j.out
#SBATCH -J VarCall_${SampleName}
hostname
date

cd ${VarDir}

echo "Sample : ${SampleSource}"
module load freebayes/1.3.4

if [ "${SampleSource}" == "WW" ]; then
	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -C 5 --pooled-continuous -b ${MapDir}/minimap_${SampleName}_PE.bam > ${SampleName}_PE.vcf
else
	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -F 0.2 -b ${MapDir}/minimap_${SampleName}_PE.bam > ${SampleName}_PE.vcf
fi

#snpEff : Wuhan-Hu-1_NC_045512.2
module load snpEff/4.3q
java -Xmx20g -jar /isg/shared/apps/snpEff/4.3q/snpEff.jar Wuhan-Hu-1_NC_045512.2 ${SampleName}_PE.vcf > ${SampleName}_PE_Anno.vcf
echo "$0 $@"

rm -r ${SampleDir}/PE
	
VARCALL`





