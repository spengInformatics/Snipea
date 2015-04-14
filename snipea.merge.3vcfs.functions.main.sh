#!/bin/bash

function splitSNVsINDELsFromSeuratFile(){
	echo -e "Splitting Seurat file up into snvs AND indels" 1>&2 ;
	local F=$1
	local HEADER=`grep -E "^#" ${F}`
	local SNVS=`grep -vE "^#" ${F} | awk -F"\t" 'length($4)<2 && length($5)<2'`
	local INDELS=`grep -vE "^#" ${F} | awk -F"\t" 'length($4)>1 || length($5)>1'`
	SEU_SNV=`basename ${F} ".vcf"`.snvs.vcf
	SEU_INDEL=`basename ${F} ".vcf"`.indels.vcf
	if [[ -e ${SEU_SNV} ]] ; then rm ${SEU_SNV} ; fi
	if [[ -e ${SEU_INDEL} ]] ; then rm ${SEU_INDEL} ; fi
	echo -e "${HEADER}\n${SNVS}" > ${SEU_SNV}
	echo -e "${HEADER}\n${INDELS}" > ${SEU_INDEL}
}

function captureDuplicatedPositions(){
	local LISTDUPPOS=""
	for F in $* ; do 
		if [[ "${LISTDUPPOS[@]}" != "" ]] ; then 
			LISTDUPPOS="`cat $F | grep -vE \"^#\" | cut -f1-2 | sort | uniq -d | awk '{if($0 != "" ) { print $1"\t"$2-1"\t"$2 } else {printf ""} }'  | sed 's/ \+//'`"; 
		else 
			LISTDUPPOS="`cat $F | grep -vE \"^#\" | cut -f1-2 | sort | uniq -d | awk '{if($0 != "" ) { print $1"\t"$2-1"\t"$2 } else {printf ""} }'  | sed 's/ \+//'` ${LISTDUPPOS[@]}" 
		fi
	done
	if [[ "${LISTDUPPOS[@]}" != "" ]] ; then echo -e "${LISTDUPPOS[@]}" ; else echo -e "" ; fi
}

function captureDuplicatedPositionsIndels(){
	rm duplicated_positions.indels.bed >/dev/null 2>&1
	for F in $* ; do cat $F | grep -vE "^#" | cut -f1-2 | sort | uniq -d | awk '{if($0 != "" ) { print $1"\t"$2-1"\t"$2 } else {printf ""} }' | sed 's/ \+//g ; s/\t\t//g' >> duplicated_positions.indels.bed  ; done
}
function captureDuplicatedPositionsSnvs(){
	rm duplicated_positions.snvs.bed > /dev/null 2>&1
	for F in $* ; do cat $F | grep -vE "^#" | cut -f1-2 | sort | uniq -d | awk '{if($0 != "" ) { print $1"\t"$2-1"\t"$2 } else {printf ""} }'| sed 's/ \+//g ; s/\t\t//g' >> duplicated_positions.snvs.bed  ; done
}

function getPrefixOutput(){
	## $1 is PREFIX, $2 is a FILE from which the prefix should be extracted if ${PREFIX} is EMPTY
	local PREFIX=$1 ; local F=$2
	if [[ ${PREFIX} == "" ]] ; then PREFIX="`echo -e "${F}" | cut -d"." -f1`" ; fi
	echo -e "${PREFIX}" ; exit
}

function checkFile(){
	local FLAG="OK";
	for F in $* ; 
	do 
		if [[ ( ! -e $F && ! -L $F ) || ${F} == ""  ]] ; 
		then
			usage 1>&2
			echo -e "\n************************\nFILE NOT FOUND << ${F} >>. \n************************\n" 1>&2;
	#change FLAG to "File Not Found(FNF)"
			FLAG="FNF"
		fi ; 
	done
	if [[ "${FLAG}" == "OK" ]] ; then echo -e "${FLAG}" >/dev/null 2>&1 ; else echo -e "${FLAG}" ; fi
}


function countEventsPerFile(){
	## $1 is the prefix_output
	## then $* represent all the VCF files we want to grab the count of events from
	local threshold=10000  ## Threshold of number of calls in vcf; if above, we send a warning
	local prefix=$1; shift
	arrWarning=()
	for F in $* ; do 
		COUNT=`grep -vE "^#" ${F} | wc -l`
		echo -e "$F  \t${COUNT}" #| tee -a ${prefix}.variants_counts.txt
		if [[ ${COUNT} -gt ${threshold} ]] ; then arrWarning=( ${arrWarning[@]}  "${F}  ---> ${COUNT}" ) ; fi
	done
	if [[ ${arrWarning[@]} != "" ]] ; 
	then 
		echo -e "One or more files contains a high number of events.\nProcess of SNIPEA might take a while...\n`for ((i=0;i<${#arrWarning[@]};i=i+1)) ; do echo -e \"${arrWarning[$i]}\" ; done`" | mail -s "WARNING from MERGER 3 VCFs - ${prefix}" ${LOGNAME}@tgen.org
	fi
}


function cleaningTempFiles(){
	##$1 is the logical value telling us about keeping Venn Input file or not?
	keepVennInputFiles=$1
	rm CUR_LOCUS* bedtools.multiintersection*  countsOfCommon* header.vcf tmp.*  >/dev/null 2>&1
	rm L1.manage.indels.bed  L2.manage.indels.bed leftover.slk.bed slk.dup.already_processed.txt dup_indel_already_processed_as_dup.txt dup_snv_already_processed_as_dup.txt >/dev/null 2>&1
	rm duplicated_positions.snvs.bed duplicated_positions.indels.bed duplicated_positions.both_snvs_and_indels.bed  >/dev/null 2>&1
	if [[ ${keepVennInputFiles} != "yes" ]]
	then
		rm R_VennDiagram*.ini slk.venn seu.venn mtc.venn seu.indels.venn slk.indels.venn  >/dev/null 2>&1
	fi
}
