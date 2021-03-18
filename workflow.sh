
OutDir="Path/to/outputDIr"
ProjectDir="Path/to/ProjectDir"
SampleName="SampleName"
SampleDir=${ProjectDir}/${SampleName}
SampleSource="SampleSourceType WW/SI"
ReadType="PE/SE/PEtaSE"
ProjectID="PID"

queue=general
QOS=general
date

logTimeStamp=`date +%d"_"%b"_"%H"hr-"%M"min-"%S"sec"`

LogDir=${OutDir}/log-${logTimeStamp}
fastqcDir
trimDir
MapDir
VarDir




if [ ! -d "${LogDir}" ]; then
    mkdir "${LogDir}"
    chmod 777 "${LogDir}"
fi

### Merging Files
cd ${SampleDir}

if [ ${ReadType} = "PE" ]
then
		jid1=`sbatch <<- MERGE | egrep -o -e "\b[0-9]+$"
		        #!/bin/bash -l
		        #SBATCH -p ${queue}
		        #SBATCH -q ${QOS}
		        #SBATCH -c 1
		        #SBATCH -o ${LogDir}/merge-%j.out
		        #SBATCH -e $debugdir/merge-%j.err
		        #SBATCH -J PE_merge_files
		        hostname
		        date

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

		        echo "$0 $@"
		MERGE`
elif [ ${ReadType} = "PEtaSE" ]; then
		jid1=`sbatch <<- MERGE | egrep -o -e "\b[0-9]+$"
		        #!/bin/bash -l
		        #SBATCH -p ${queue}
		        #SBATCH -q ${QOS}
		        #SBATCH -c 1
		        #SBATCH -o ${LogDir}/merge-%j.out
		        #SBATCH -e $debugdir/merge-%j.err
		        #SBATCH -J PEtaSE
		        hostname
		        date

		        R1files=`ls *R1*.fastq.gz 2>/dev/null | wc -l`
		        R2files=`ls *R2*.fastq.gz 2>/dev/null | wc -l`


		        	echo -e "There is ${R1files} input file for R1 and ${R2files} for R2"
		        	echo -e "Files will be merged to a single file"
		        	zcat `ls *.fastq.gz ` >> ${SampleName}_PEtaSE.fastq


		        echo "$0 $@"
		MERGE`
else [ ${ReadType} = "SE" ]; then
		jid1=`sbatch <<- MERGE | egrep -o -e "\b[0-9]+$"
		        #!/bin/bash -l
		        #SBATCH -p ${queue}
		        #SBATCH -q ${QOS}
		        #SBATCH -c 1
		        #SBATCH -o ${LogDir}/merge-%j.out
		        #SBATCH -e $debugdir/merge-%j.err
		        #SBATCH -J SE_merge
		        hostname
		        date

		        Readfiles=`ls *.fastq.gz 2>/dev/null | wc -l`
		        	echo -e "There is ${Readfiles} input fastq.gz"
		        	echo -e "Files will be merged to a single file"
		        	zcat `ls *.fastq.gz ` >> ${SampleName}_SE.fastq


		        echo "$0 $@"
		MERGE`
fi

dependmerge="afterok:$jid1"


















