#!/bin/bash
#USAGE: Rank all SNV based on weighted score (related to Callers agreement, SNV effect, callers quality score and public database status)

function rankSNVs(){
	local snpeffFile=$1
    rm -f tmp.final_RANKED.snvs.tsv
	rm -f tmp.final_RANKED_SORTED.snvs.tsv
	while read -r LINE ; do 
        weightedScore=0
#	    COUNT_NUM=`echo ${LINE}|cut -d" " -f8|cut -d";" -f1|cut -d"=" -f2`
		COUNT_NUM_TAG=`echo ${LINE}|cut -d" " -f8|awk -F";" '{for(i=1;i<=NF;i++){ if($i ~ /^CALLERS_COUNT=/) { print i ; break }}}'`
	    COUNT_NUM=`echo ${LINE}|cut -d" " -f8|cut -d";" -f${COUNT_NUM_TAG}|cut -d"=" -f2`
		
		CALLERS_QS=`echo ${LINE}|cut -d" " -f6`
		DATABASE_STATUS=`echo ${LINE}|grep COSMIC`
		SNPEFF_STATUS=`echo ${LINE}|egrep -h "HIGH|NON_SYNONYMOUS_CODING"`

	    case ${COUNT_NUM} in
			1) weightedScore=$((${weightedScore}+100));;
			2) weightedScore=$((${weightedScore}+200));;
			3) weightedScore=$((${weightedScore}+300));;
			*) exit ;;
		esac
		
		if [[ ${DATABASE_STATUS} != "" ]]  
	       then 
		 weightedScore=$((${weightedScore}+40))
		fi
		
		if [[ ${SNPEFF_STATUS} != "" ]]  
	       then 
		if [[ ${COUNT_NUM} -eq 1 ]];  then weightedScore=$((${weightedScore}+50));else weightedScore=$((${weightedScore}+100));fi
		fi
		
		weightedScore=`echo "${weightedScore}+${CALLERS_QS}"|bc`
        weightedScore=`printf "%.2f" $(echo "${weightedScore}/5"|bc -l)`
		echo ${LINE}|awk -v weightedScore="${weightedScore}" '{for(i=1;i<=NF;i++) {if ( i==7 ) {printf weightedScore"\t";i++} else if ( i==NF ) {printf $i"\n";break} {printf $i"\t";}};}' >>  tmp.final_RANKED.snvs.tsv
	done < ${snpeffFile}
    cat tmp.final_RANKED.snvs.tsv|sort -r -k7,7 > tmp.final_RANKED_SORTED.snvs.tsv
}

AnnotatedFile=$1
cd $2
PREFIX=`basename ${AnnotatedFile}|cut -d. -f1`
cat ${AnnotatedFile}|grep "^#" > tmp.new_header.vcf
inputFile=tmp.${PREFIX}.beforeRank.vcf
cat ${AnnotatedFile}|grep -vE "^#" > ${inputFile}
rankSNVs "${inputFile}"
cat tmp.new_header.vcf  tmp.final_RANKED_SORTED.snvs.tsv | sed '/^$/d' >  "${PREFIX}.final.ranked.vcf"
exit
