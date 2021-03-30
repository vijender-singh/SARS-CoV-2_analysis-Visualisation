#!/bin/bash

configFilePath=${1}

if [ -z ${configFilePath} ];then
	source ./config_file.txt
else
	source ${configFilePath}/config_file.txt
fi

# Variable check

#mkdir -p ${OutDir}/${ProjectID}

logTimeStamp=`date +%d"_"%b"_"%H"hr-"%M"min-"%S"sec"`
ProjectLogDir=${OutDir}/log-${ProjectID}
OUTDIR=${OutDir}/${ProjectID}_${AnalysisType}
LogDir=${ProjectLogDir}/log-${logTimeStamp}_${AnalysisType}

if [ ! -d "${LogDir}" ]; then
    mkdir -p "${LogDir}"
    chmod 777 "${LogDir}"
fi


if [ -z ${MultiTypeAnalysis}];then
		lensampleNames=${#SampleNames[@]}
		lenSampleSources=${#SampleSources[@]}
		if [ ${lensampleNames} = ${lenSampleSources} ]; then
				for indx in `seq 1 ${lensampleNames}`
				do
				sbatch -o ${LogDir}/mainscript_${SampleNames[ `expr ${indx} - 1` ]}_log.out ./SLURM/workflow_${AnalysisType}.sh ${SampleNames[ `expr ${indx} - 1` ]} ${SampleSources[ `expr ${indx} - 1` ]} ${DataDir} ${OUTDIR} ${ProjectID} ${ResourceDir} ${logTimeStamp} ${LogDir} ${ProjectLogDir}
				done
		 else
		 		echo -e " ERROR : The Number of Samples  is not equal to number of Sources.  Make sure that each sample has its corresponding source listed."
		 fi

elif [ ${#MultiTypeAnalysis[@]} -gt 0 ]; then
 		lenMultiAnalysis=${#MultiTypeAnalysis[@]}
		lensampleNames=${#SampleNames[@]}
		lenSampleSources=${#SampleSources[@]} 
		if [ ${lensampleNames} = ${lenSampleSources} ] && [ ${lensampleNames} = ${lenMultiAnalysis} ] ; then
				for indx in `seq 1 ${lensampleNames}`
				do
				sbatch -o ${LogDir}/mainscript_${SampleNames[ `expr ${indx} - 1` ]}_log.out ./SLURM/workflow_${AnalysisType}.sh ${SampleNames[ `expr ${indx} - 1` ]} ${SampleSources[ `expr ${indx} - 1` ]} ${DataDir} ${OUTDIR} ${ProjectID} ${ResourceDir} ${logTimeStamp} ${LogDir} ${ProjectLogDir}
				done
		 else
		 		echo -e and " ERROR : listed The listed Number of Samples,  Sources and analysis type donot match.  Make sure that each sample has its corresponding source and its analysis type is listed."
		 fi

else
	
	echo -e "ERROR : Config file has issues, please check it again and make appropriate corrections."

fi

         
