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
TEMP_INDEL_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
SEUBED=$3	## mandatory BED file - BED format expected
SLKBED=$4	## mandatory BED file - BED format expected

function processSEUonly(){
	old_IFS=$IFS		# save the field separator
	IFS=$'\t'
	new_IFS=${IFS}
	local CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
	local TEMP_INDEL_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
	local SEUBED=$3	## mandatory BED file - BED format expected
	CALLERS_COUNT=1; FLAG_TOOL="SEURAT";
	while ${IFS} read LINE ; do
		MCHR=`echo -e "${LINE}" | cut -f1`
		MPOS=`echo -e "${LINE}" | cut -f3`
		MID="."
		MREF=`echo -e "${LINE}" | cut -f5`
		MALT=`echo -e "${LINE}" | cut -f6`
		MQUAL=`echo -e "${LINE}" | cut -f7`
		MFILT=`echo -e "${LINE}" | cut -f8`
		CALLERS_COUNT=1; FLAG_TOOL="SEURAT";
		local L1=(`echo -e "${LINE}" | cut -f1,3-`)
		MINFO=$(concatInfoFields "${L1[7]}")
		MFORMAT=$(indels_getFORMAT "${L1[7]}" "" "" "" "" "${MREF}" "${MALT}" "SEURAT_DUPLICATE")
		echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO}${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
	done < ${SEUBED}
	IFS=${old_IFS}
}


function processSLKonly(){
	old_IFS=$IFS		# save the field separator
	IFS=$'\t' 
	new_IFS=${IFS}
	local CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
	local TEMP_INDEL_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
	local SLKBED="$3"	## mandatory BED file - BED format expected
	CALLERS_COUNT=1; FLAG_TOOL="STRELKA";
	if [[ -s ${SLKBED} ]]
	then 
	while ${IFS} read LINE ; do
		MCHR=`echo -e "${LINE}" | cut -f1`
		MPOS=`echo -e "${LINE}" | cut -f3`
		MID="."
		MREF=`echo -e "${LINE}" | cut -f5`
		MALT=`echo -e "${LINE}" | cut -f6`
		MQUAL=`echo -e "${LINE}" | cut -f7`
		MFILT=`echo -e "${LINE}" | cut -f8`
		CALLERS_COUNT=1; FLAG_TOOL="STRELKA";
		local L1=(`echo -e "${LINE}" | cut -f1,3-`)
		MINFO=$(concatInfoFields "${L1[7]}")
		MFORMAT=$(indels_getFORMAT "${L1[7]}" "" "" "" "" "${MREF}" "${MALT}" "SLK")
		echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
	done < ${SLKBED}
	fi
	IFS=${old_IFS}
}



function processSEUandSLK(){
	old_IFS=$IFS		# save the field separator
	IFS=$'\t' 
	new_IFS=${IFS}
	rm slk.dup.already_processed.txt > /dev/null 2>&1
	local CUR_LOCUS_FILE="$1"	## mandatory BED file - BED format expected
	local TEMP_INDEL_OUTFILE="$2"	### tabulated text file following vcf specifiction (but without the vcf header; will be used later for concatenation)
	local SEUBED=$3	## mandatory BED file - BED format expected
	local SLKBED=$4	## mandatory BED file - BED format expected

	TMPL1=L1.manage.indels.bed 	#L1 to manage SEU data
	TMPL2=L2.manage.indels.bed	#L2 to manage SLK data
	bedtools intersect -a ${SEUBED} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL1} &
	bedtools intersect -a ${SLKBED} -b ${CUR_LOCUS_FILE} -f 1 > ${TMPL2} &
	wait

	while ${IFS} read LINE ; do
		L1=() ; L2=() ; 
		CALLERS_COUNT=1; 
		FLAG_TOOL="SEURAT";
		PREVALENCE="SEU"
		FLAG_COMMON_POS="NO"
		MCHR=`echo -e "${LINE}" | cut -f1`
		MPOS=`echo -e "${LINE}" | cut -f3`
		MID="."
		MREF=`echo -e "${LINE}" | cut -f5` ; local MREF_SEU=`echo -e "${LINE}" | cut -f5`
		MALT=`echo -e "${LINE}" | cut -f6` ; local MALT_SEU=`echo -e "${LINE}" | cut -f6`
		REBUILTPOS_SEU="${MCHR}${MPOS}${MREF}"
		## now we need to check if the same position with the same ref is present in SLK as well"
		while ${IFS} read LINESLK ; do
			MCHRSLK=`echo -e "${LINESLK}" | cut -f1`
			MPOSSLK=`echo -e "${LINESLK}" | cut -f3`
			MIDSLK="."
			MREFSLK=`echo -e "${LINESLK}" | cut -f5`
			MALTSLK=`echo -e "${LINESLK}" | cut -f6`
			REBUILTPOS_SLK="${MCHRSLK}${MPOSSLK}${MREFSLK}"
			CC=0 ### count the common position between SEU and SLK for a specific position
			if [[ "${REBUILTPOS_SEU}" == "${REBUILTPOS_SLK}" ]] ; 
			then
				CALLERS_COUNT=2; FLAG_TOOL="SEURAT;STRELKA"; PREVALENCE="SLK"
				CC=$((CC+1)) ; 
				FLAG_COMMON_POS="YES" ;
		 		L1=(`echo -e "${LINE}" | cut -f1,3-`)
				L2=(`echo -e "${LINESLK}" | cut -f1,3-`)
				MREF=$(indels_getREF "${L1[3]}" "${L2[3]}")
				MALT=$(indels_getALT "${L1[4]}" "${L2[4]}")
				if [[ "${PREVALENCE}" == "SEU" ]] ; then MQUAL=${L1[5]} ; else MQUAL="." ; fi
				MFILT="."
				MINFO=$(concatInfoFields "${L1[7]}" "${L2[7]}")
				MFORMAT=$(indels_getFORMAT "${L1[7]}" "${L2[7]}" "${L2[8]}" "${L2[9]}" "${L2[10]}" "${MREF}" "${MALT}" "${PREVALENCE}")
				echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
				grep -E "^${MCHR}" ${TMPL2} | grep -E "${MPOS}" | awk -v oREF=${MREF} -v oALT=${MALT} '$5==oREF && $6==oALT ' >> slk.dup.already_processed.txt
				IFS=${new_IFS}
			else
				## According to JK's feedback: This time if we have a call of an indel from both seurat and strelka, the prevalence goes to Strelka even though the REF is different from the given tools, we keep STRELKA, but I added the "REF" and ALT call from SEURAT as well.
				## Actually if the REF (or ALT) is different in Strelka and Seurat, we cretate two lines, so this section shold never happens
				CALLERS_COUNT=2; FLAG_TOOL="SEURAT;STRELKA"; PREVALENCE="SLK"
				FLAG_COMMON_POS="YES" ;
				L1=(`echo -e "${LINE}" | cut -f1,3-`)
				L2=(`echo -e "${LINESLK}" | cut -f1,3-`)
				MREF=$(indels_getREF "" "${L2[3]}")
				MALT=$(indels_getALT "" "${L2[4]}")
				MQUAL="."
				MFILT="." # "${L2[6]}"  otherwise
				MINFO=$(concatInfoFields "${L1[7]}" "${L2[7]}")
				MFORMAT=$(indels_getFORMAT "" "${L2[7]}" "${L2[8]}" "${L2[9]}" "${L2[10]}" "${MREF}" "${MALT}" "${PREVALENCE}")
				echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};SEURAT_REF=${MREF_SEU};SEURAT_ALT=${MALT_SEU};${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
				grep -E "^${MCHR}" ${TMPL2} | grep -E "${MPOS}" | awk -v oREF=${MREF} -v oALT=${MALT} '$5==oREF && $6==oALT ' >> slk.dup.already_processed.txt
				IFS=${new_IFS}
			fi
		done < ${TMPL2}

		if [[ ${FLAG_COMMON_POS} == "YES" && -s slk.dup.already_processed.txt ]]
		then 
			bedtools intersect -a ${TMPL2} -b slk.dup.already_processed.txt -v > leftover.slk.bed 
			if [[ `wc -l leftover.slk.bed | awk '{print $1}'` -ge 1 ]] ; then processSLKonly "${CUR_LOCUS_FILE}" "${TEMP_INDEL_OUTFILE}" "leftover.slk.bed" ; fi
			continue ;
		else
			continue
		fi
		MQUAL=`echo -e "${LINE}" | cut -f7`
		MFILT=`echo -e "${LINE}" | cut -f8`
		L1=(`echo -e "${LINE}" | cut -f1,3-`)
		MINFO=$(concatInfoFields "${L1[7]}")
		MFORMAT=$(indels_getFORMAT "${L1[7]}" "" "" "" "" "${MREF}" "${MALT}" "${PREVALENCE}")
		echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
		## then we check if Strelka has still lines in the file ${SLKBED} ; if so we process it with   processSLKonly
		if [[ -s slk.dup.already_processed.txt ]] ; then  bedtools intersect -a ${TMPL2} -b slk.dup.already_processed.txt -v > leftover.slk.bed ; fi
		if [[ `wc -l leftover.slk.bed | awk '{print $1}'` -ge 1 ]] ; then processSLKonly "${CUR_LOCUS_FILE}" "${TEMP_INDEL_OUTFILE}" "leftover.slk.bed" ; fi

	done < ${TMPL1}
	IFS=${old_IFS}

}


###################
###### MAIN #######
###################
nrowSEU=`cat "${SEUBED}" | grep -vE "^#" | wc -l`
nrowSLK=`cat "${SLKBED}" | grep -vE "^#" | wc -l`

if [[ ${nrowSEU} -ge 1 && ${nrowSLK} -eq 0 ]] ; then processSEUonly "${CUR_LOCUS_FILE}" "${TEMP_INDEL_OUTFILE}" "${SEUBED}" ; fi
if [[ ${nrowSLK} -ge 1 && ${nrowSEU} -eq 0 ]] ; then processSLKonly "${CUR_LOCUS_FILE}" "${TEMP_INDEL_OUTFILE}" "${SLKBED}" ; fi
if [[ ${nrowSEU} -ge 1 && ${nrowSLK} -ge 1 ]] ; then processSEUandSLK "${CUR_LOCUS_FILE}" "${TEMP_INDEL_OUTFILE}" "${SEUBED}" "${SLKBED}" ; fi
exit

