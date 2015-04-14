#!/bin/bash
## THIS SCRIPT ONLY MANAGES ONE LOCUS.
## example of CUR_LOCUS_FILE: 11	59836985	59836986
## NOTE: if this script is run it is because AT LEAST one POSITION has been found twice in the VCF file. So we expect that at least SEURAT or STRELKA files have duplicate;
## therefore the first loop dealing with seurat or the second loop dealing with strelak should run at least ONCE.
#bash ${dir_script}/snipea.merge.3vcfs.source_functions.sh ${dir_script}
dir_script=`dirname $0`
source ${dir_script}/snipea.merge.3vcfs.functions.snvs.sh
source ${dir_script}/snipea.functions.static.sh
source ${dir_script}/snipea.merge.3vcfs.functions.indels.sh

CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
TEMP_SNV_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
SEUBED=$3	## mandatory BED file - BED format expected
SLKBED=$4	## mandatory BED file - BED format expected
MTCBED=$5	## mandatory BED file - BED format expected

function processOneCallerOnly(){
	old_IFS=$IFS		# save the field separator
	IFS=$'\t'
	new_IFS=${IFS}
	local CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
	local TEMP_SNV_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
	local SEUBED=$3	## mandatory BED file - BED format expected
	local FLAG_TOOL=$4
	CALLERS_COUNT=1;
	while ${IFS} read LINE ; do
		MCHR=`echo -e "${LINE}" | cut -f1`
		MPOS=`echo -e "${LINE}" | cut -f3`
		MID="."
		MREF=`echo -e "${LINE}" | cut -f5`
		MALT=`echo -e "${LINE}" | cut -f6`
		MQUAL=`echo -e "${LINE}" | cut -f7`
		MFILT=`echo -e "${LINE}" | cut -f8`
		local L1=(`echo -e "${LINE}" | cut -f1,3-`)
		MINFO=$(concatInfoFields "${L1[7]}")
		MFORMAT=$(indels_getFORMAT "${L1[7]}" "" "" "" "" "${MREF}" "${MALT}" "${FLAG_TOOL}_DUPLICATE")
		echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_SNV_OUTFILE}
	done < ${SEUBED}
	IFS=${old_IFS}
}

function processTwoCallers(){
	old_IFS=$IFS		# save the field separator
	IFS=$'\t' 
	new_IFS=${IFS}
	rm slk.dup.already_processed.txt > /dev/null 1>&2
	local CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
	local TEMP_SNV_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
	local TOOL1=$3	## mandatory BED file - BED format expected
	local TOOL2=$4	## mandatory BED file - BED format expected
	local FLAG_TOOL=$5
	local PREVALENCE=$6

	TMPL1=L1.manage.indels.bed 	#L1 to manage SEU data
	TMPL2=L2.manage.indels.bed	#L2 to manage SLK data
	bedtools intersect -a ${TOOL1} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL1}
	bedtools intersect -a ${TOOL2} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL2}

	L1=(`cat "${TMPL1}" | cut -f1,3-`)
	L2=(`cat "${TMPL2}" | cut -f1,3-`)
	MCHR=${L1[0]}
	MPOS=${L1[1]}
	MID="."
	MREF=$(getREF "${L1[3]}" "${L2[3]}" "${L3[3]}")
	MALT=$(getALT "${L1[4]}" "${L2[4]}" "${L3[4]}")
	MQUAL=$(getQUAL "${L1[5]}" "${L2[5]}" "${L3[5]}")
	MFILT="."
	MINFO=$(concatInfoFields "${L1[7]}" "${L2[7]}" "${L3[7]}")
	MFORMAT=$(getFORMAT "${L1[7]}" "${L2[7]}" "${L2[8]}" "${L2[9]}" "${L2[10]}" "${L3[7]}" "${L3[8]}" "${L3[9]}" "${L3[10]}" "${MREF}" "${MALT}" "${PREVALENCE}")
	echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_SNV_OUTFILE}

}

function processThreeCallers(){
	old_IFS=$IFS		# save the field separator
	IFS=$'\t' 
	new_IFS=${IFS}
	rm slk.dup.already_processed.txt > /dev/null 1>&2
	local CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
	local TEMP_SNV_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
	local TOOL1=$3	## mandatory BED file - BED format expected
	local TOOL2=$4	## mandatory BED file - BED format expected
	local TOOL3=$5	## mandatory BED file - BED format expected
	local FLAG_TOOL=$6
	local PREVALENCE=$7

	TMPL1=L1.manage.indels.bed 	#L1 to manage SEU data
	TMPL2=L2.manage.indels.bed	#L2 to manage SLK data
	TMPL3=L3.manage.indels.bed	#L2 to manage MTC data
	bedtools intersect -a ${TOOL1} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL1}
	bedtools intersect -a ${TOOL2} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL2}
	bedtools intersect -a ${TOOL3} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL3}

	L1=(`cat "${TMPL1}" | cut -f1,3-`)
	L2=(`cat "${TMPL2}" | cut -f1,3-`)
	L3=(`cat "${TMPL3}" | cut -f1,3-`)
	MCHR=${L1[0]}
	MPOS=${L1[1]}
	MID="."
	MREF=$(getREF "${L1[3]}" "${L2[3]}" "${L3[3]}")
	MALT=$(getALT "${L1[4]}" "${L2[4]}" "${L3[4]}")
	MQUAL=$(getQUAL "${L1[5]}" "${L2[5]}" "${L3[5]}")
	MFILT="."
	MINFO=$(concatInfoFields "${L1[7]}" "${L2[7]}" "${L3[7]}")
	MFORMAT=$(getFORMAT "${L1[7]}" "${L2[7]}" "${L2[8]}" "${L2[9]}" "${L2[10]}" "${L3[7]}" "${L3[8]}" "${L3[9]}" "${L3[10]}" "${MREF}" "${MALT}" "${PREVALENCE}")
	echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_SNV_OUTFILE}
}

###################
###### MAIN #######
###################
nrowSEU=`cat "${SEUBED}" | grep -vE "^#" | wc -l`
nrowSLK=`cat "${SLKBED}" | grep -vE "^#" | wc -l`
nrowMTC=`cat "${MTCBED}" | grep -vE "^#" | wc -l`

if [[ ${nrowSEU} -ge 1 && ${nrowSLK} -eq 0 && ${nrowMTC} -eq 0 ]] ; then processOneCallerOnly "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${SEUBED}" "SEURAT"; fi
if [[ ${nrowSLK} -ge 1 && ${nrowSEU} -eq 0 && ${nrowMTC} -eq 0 ]] ; then processOneCallerOnly "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${SLKBED}" "STRELKA" ; fi
if [[ ${nrowMTC} -ge 1 && ${nrowSEU} -eq 0 && ${nrowSLK} -eq 0 ]] ; then processOneCallerOnly "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${MTCBED}" "MUTECT" ; fi
if [[ ${nrowSEU} -ge 1 && ${nrowSLK} -ge 1 ]] ; then processTwoCallers "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${SEUBED}" "${SLKBED}" "SEURAT;STRELKA" "SEU"; fi
if [[ ${nrowSEU} -ge 1 && ${nrowMTC} -ge 1 ]] ; then processTwoCallers "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${SEUBED}" "${SLKBED}" "SEURAT;MUTECT" "SEU"; fi
if [[ ${nrowSLK} -ge 1 && ${nrowMTC} -ge 1 ]] ; then processTwoCallers "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${SLKBED}" "${MTCBED}" "STRELKA;MUTECT" "SLK"; fi
if [[ ${nrowSEU} -ge 1 && ${nrowSLK} -ge 1 && ${nrowMTC} -ge 1 ]] ; then processThreeCallers "${CUR_LOCUS_FILE}" "${TEMP_SNV_OUTFILE}" "${SEUBED}" "${SLKBED}" ${MUTECT} "SEURAT;STRELKA;MUTECT" "SEU"; fi
exit

