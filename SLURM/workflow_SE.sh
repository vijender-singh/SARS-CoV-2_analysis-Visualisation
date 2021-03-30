#!/bin/bash
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
mkdir -p SE

rm *.fastq

Readfiles="$(ls -1 *.fastq.gz | wc -l)"

if [ -f ./SE/${SampleName}_SE.fastq ]; then
        rm ./SE/${SampleName}_SE.fastq
fi

zcat `ls *R*.fastq.gz` >> SE/${SampleName}_SE.fastq


reads=$(wc -l < ./SE/${SampleName}_SE.fastq)

total_reads=$(($reads / 4))

if [ ${total_reads} -eq 0 ]; then
	echo -e "ABORT ABORT ABORT !!!!....   No reads detected in merged file.   Please check !!!"
	exit 1
fi



echo -e "Total Reads : ${reads}"



jid2=`sbatch --parsable <<- TRIMMOTRIM
#!/bin/bash -l
#SBATCH -p ${queue}
#SBATCH -q ${qos}
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -o ${LogDir}/${SampleName}-QC_TRIM_QC-SE_%j.out
hostname
date
cd ${SampleDir}
module load fastqc/0.11.7

fastqc -o ${RawfastqcDir} --threads 4  ./SE/${SampleName}_SE.fastq
module rm fastqc/0.11.7
module load Trimmomatic/0.39

java -jar /isg/shared/apps/Trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 4 \
./SE/${SampleName}_SE.fastq \
${trimDir}/trim_${SampleName}_SE.fastq \
ILLUMINACLIP:${resourceDir}/TruSeq3-SE.fa:2:30:10 \
SLIDINGWINDOW:4:25 \
MINLEN:45

module rm Trimmomatic/0.39
module load fastqc/0.11.7

fastqc -o ${TrimfastqcDir} --threads 4  ${trimDir}/trim_${SampleName}_SE.fastq

module rm fastqc/0.11.7
		
echo "$0 $@"
TRIMMOTRIM`



depend2="afterok:${jid2}"
echo depend2

jid3=`sbatch --parsable <<- MAPPING
#!/bin/bash -l
#SBATCH -p ${queue}
#SBATCH -q ${qos}
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -d ${depend2}
#SBATCH -o ${LogDir}/${SampleName}-map-SE_%j.out
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

samtools depth -d 0 -Q 30 -q 0 -aa minimap_${SampleName}_SE.bam > ${VarDir}/Depth_minimap_${SampleName}_SE.txt

bedtools genomecov -bg -ibam minimap_${SampleName}_SE.bam >> ${VarDir}/${SampleName}_SE.bdg

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
#SBATCH -o ${LogDir}/${SampleName}-VarCall-SE_%j.out
#SBATCH -J VarCall_${SampleName}
hostname
date

cd ${VarDir}

echo "Sample : ${SampleSource}"
module load freebayes/1.3.4

if [ "${sampleSource}" = "WW" ]; then
	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -C 5 --pooled-continuous ${MapDir}/minimap_${SampleName}_SE.bam > ${SampleName}_SE.vc
else
	freebayes -f ${resourceDir}/Wuhan-Hu-1_NC_045512.2.fasta -m 30 -q 20 -F 0.2 ${MapDir}/minimap_${SampleName}_SE.bam > ${SampleName}_SE.vcf
fi

#snpEff : Wuhan-Hu-1_NC_045512.2
module load snpEff/4.3q

java -Xmx20g -jar /isg/shared/apps/snpEff/4.3q/snpEff.jar Wuhan-Hu-1_NC_045512.2 ${SampleName}_SE.vcf > ${SampleName}_SE_Anno.vcf

rm -r ${SampleDir}/SE

echo "$0 $@"
VARCALL`























