#!/bin/bash

source config_file.txt

# Variable check


if [ -z ${MultiTypeAnalysis}];then
		lensampleNames=${#SampleNames[@]}
		lenSampleSources=${#SampleSources[@]}
		if [ ${lensampleNames} = ${lenSampleSources} ]; then
				for indx in `seq 1 ${lensampleNames}`
				do
				sh workflow_${AnalysisType} ${SampleNames[ `expr ${indx} - 1` ]} ${SampleSources[ `expr ${indx} - 1` ]} ${ProjectDir} ${OutDir} ${ProjectID}
		 else
		 		echo -e " ERROR : The Number of Samples  is not equal to number of Sources.  Make sure that each sample has its corresponding source listed."
		 fi

elif [ ${#MultiTypeAnalysis[@]} -gt 0 ]; then
 		lenMultiAnalysis=${#MultiTypeAnalysis[@]}
		lensampleNames=${#SampleNames[@]}
		lenSampleSources=${#SampleSources[@]} 
		if [ ${lensampleNames} = ${lenSampleSources} ] && ; then [ ${lensampleNames} = ${lenMultiAnalysis} ] 
				for indx in `seq 1 ${lensampleNames}`
				do
				sh workflow_${AnalysisType} ${SampleNames[ `expr ${indx} - 1` ]} ${SampleSources[ `expr ${indx} - 1` ]} ${ProjectDir} ${OutDir} ${ProjectID}
		 else
		 		echo -e and " ERROR : listed The listed Number of Samples,  Sources and analysis type donot match.  Make sure that each sample has its corresponding source and its analysis type is listed."
		 fi

else
	
	echo -e "ERROR : Config file has issues, please check it again and make appropriate corrections."

fi

         
