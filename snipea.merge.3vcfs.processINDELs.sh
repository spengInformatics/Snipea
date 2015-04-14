#!/bin/bash
##USAGE: Continue Snipea and start dealing with INDELS
dir_script=`dirname $0`
source ${dir_script}/snipea.merge.3vcfs.functions.snvs.sh
source ${dir_script}/snipea.functions.static.sh
source ${dir_script}/snipea.merge.3vcfs.functions.indels.sh

function getMergedINDELs(){
	#Order will always be: STRELKA,SEURAT
	local SEU="$1"
	local SLK="$2"
	local DUPLICATES_FILE=$3
	local TEMP_INDEL_OUTFILE="tmp.merged.indels.tsv"
	local CUR_LOCUS_FILE="CUR_LOCUS.INDELS.bed"
	local DUP_ALREADY_DONE_FILE="dup_indel_already_processed_as_dup.txt" ; >${DUP_ALREADY_DONE_FILE}

	grep -vE "^#"  ${SEU} | awk 'BEGIN {FS=OFS="\t"} {if($1 !~ /^#/ ) { $2=$2-1"\t"$2 ; print } else {print } }' | bedtools sort -i - > tmp.${SEU}.bed &
	grep -vE "^#"  ${SLK} | awk 'BEGIN {FS=OFS="\t"} {if($1 !~ /^#/ ) { $2=$2-1"\t"$2 ; print } else {print } }' | bedtools sort -i - > tmp.${SLK}.bed &
	wait
	bedtools multiinter -i tmp.${SEU}.bed tmp.${SLK}.bed | cut -f5 | sort | uniq -c  > countsOfCommonIndelsLoci.perTools.bedtools.txt &
	bedtools multiinter -i tmp.${SEU}.bed tmp.${SLK}.bed | awk 'BEGIN{ FS=OFS="\t" } { if ($3-$2 <2) {print $0 } else { for(c=0;c<$3-$2;c++){ printf "%s\t%s\t%s",$1,$2+c,$2+c+1  ; for(j=4;j<=NF;j++) { printf "\t%s",$j } ; printf "\n" } } } ' > bedtools.multiintersection.2vcfs.indels.seu.slk.bed &
	wait

	old_IFS=$IFS      # save the field separator
	IFS=$'\t'
	new_IFS=${IFS}
	while ${IFS} read LINE ; do 
		echo -e "$LINE" > ${CUR_LOCUS_FILE} ;
		MCHR=`echo -e "${LINE}" | cut -f1`
		MPOS=`echo -e "${LINE}" | cut -f3`
		MID="."
		PRESENT_IN="`echo -e "$LINE" | cut -f 5 | sed "s/,/\t/g"`" ; ## \t : because we changed the IFS  ## With PRESENT_IN we grab the appropriate column from the multi-intersection bedfile
		LOCUS="`echo -e "$LINE" | cut -f1-3`" ; 
		if [[ `bedtools intersect -a ${CUR_LOCUS_FILE} -b ${DUP_ALREADY_DONE_FILE} -f 1 2>/dev/null | wc -l  ` -eq 1  ]] ; then continue ; fi
		if [[ -s ${DUPLICATES_FILE} ]] ## test if file is not zero
		then
			## here we will manage the duplicated positions either they are in seurat, strelka or mutect
			if [[ `bedtools intersect -a ${CUR_LOCUS_FILE} -b ${DUPLICATES_FILE} -f 1 | wc -l ` -eq 1  ]] ;
			then
				bedtools intersect -a tmp.${SEU}.bed -b ${CUR_LOCUS_FILE} -f 1 > tmp.dupIndel.${SEU}.bed &
				bedtools intersect -a tmp.${SLK}.bed -b ${CUR_LOCUS_FILE} -f 1 > tmp.dupIndel.${SLK}.bed &
				wait
				bash ${dir_script}/snipea.merge.vcfs.indels.manage_duplicates.sh "${CUR_LOCUS_FILE}" "${TEMP_INDEL_OUTFILE}" "tmp.dupIndel.${SEU}.bed" "tmp.dupIndel.${SLK}.bed"
				bedtools intersect  -a ${DUPLICATES_FILE} -b ${CUR_LOCUS_FILE} >> ${DUP_ALREADY_DONE_FILE}  &
				bedtools intersect  -a ${DUPLICATES_FILE} -b ${CUR_LOCUS_FILE} -v > TT ; mv TT ${DUPLICATES_FILE}  &
				wait
				continue
			fi ; ## if YES we manage it here 
		fi
		L1=() ; L2=() ; CALLERS_COUNT=0; FLAG_TOOL="";
		for I in ${PRESENT_IN}
		do
		# get the LINE of Interest for each tool and add the LINE to an ARRAY; ## here the "cut -f1,3-" remove the extra column we have in the bed format to make it vcf format
		case ${I} in
			1) L1=(`bedtools intersect -a tmp.${SEU}.bed -b ${CUR_LOCUS_FILE} -f 1 | cut -f1,3-`) ; CALLERS_COUNT=$((${CALLERS_COUNT}+1)) ; FLAG_TOOL=$(addToolnameAsFlag "${FLAG_TOOL[@]}" "SEURAT")  ;;
			2) L2=(`bedtools intersect -a tmp.${SLK}.bed -b ${CUR_LOCUS_FILE} -f 1 | cut -f1,3-`) ; CALLERS_COUNT=$((${CALLERS_COUNT}+1)) ; FLAG_TOOL=$(addToolnameAsFlag "${FLAG_TOOL[@]}" "STRELKA")  ;;
			*) echo -e "ERROR in PRESENT_IN loop for INDELS; Aborting merge" 1>&2 ; exit ;;
		esac
		if [[ ${CALLERS_COUNT} -eq 2 ]]; then CALLERS_COUNT=3; fi
		done
		## Since we believe that Strelka Called INDELS better according to researchers, we Set the PREVALENCE to STRELKA for ALL INDEL cases ; Therefore, the case above is obsolete
		if [[ "${L2[@]}" != "" ]] ; then PREVALENCE="SLK" ; elif [[ "${L1[@]}" != "" ]] ; then PREVALENCE="SEU" ; else echo -e "" ; continue ; fi
		
		MREF=$(indels_getREF "${L1[3]}" "${L2[3]}") ### what if REF is represented by a deletion??? such as: TGGG  and ALT is T ; I have received only vcf containing SOMATIC CAlls
		MALT=$(indels_getALT "${L1[4]}" "${L2[4]}")
		MQUAL=$(getQUAL "${L1[5]}" "${L2[7]}")
		MFILT="."
		MINFO="$(concatInfoFields "${L1[7]}" "${L2[7]}")" ; #MINFO="`echo -e ${MINFO[@]} | sed 's/;$// ; s/;;/;/g'`"
		MFORMAT=$(indels_getFORMAT "${L1[7]}" "${L2[7]}" "${L2[8]}" "${L2[9]}" "${L2[10]}" "${MREF}" "${MALT}" "${PREVALENCE}")

		if [[ ! ${MREF} =~ "," ]] ; 
		then 
			echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO}${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
			## There is no ";" between ${MINFO}${MFORMAT} because of how we handle MFORMAT on the indels_getFORMAT function; This function also deals with MIF and return more than only FORMAT fields
		else
			echo -e "WARNING: this POSITION is AMBIGUOUS: ${MCHR}:${MPOS}" 1>&2 ; 
			# splitting the POSITION into two lines:
			local ARR_REF=(`echo ${MREF} | sed 's/,/\t/g'`) ; CREF=${#ARR_REF[@]} ; #echo -e "ARR_REF == ${ARR_REF[@]} \t CREF= ${CREF} ; \t ${MREF}" 1>&2 ;
			for (( L=1; L<=${CREF}; L++ ))
			do
				local newMREF="$(echo ${MREF} | cut -d',' -f ${L})"
				local newMALT=`echo -e ${MALT} | cut -d',' -f ${L} `
				echo -e "${MCHR}\t${MPOS}\t${MID}\t`echo -e \"${MREF}\" | cut -d',' -f${L}`\t`echo -e \"${MALT}\" | cut -d',' -f${L}`\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO}${MFORMAT}" >> ${TEMP_INDEL_OUTFILE}
				## There is no ";" between ${MINFO}${MFORMAT} because of how we handle MFORMAT on the indels_getFORMAT function; This function also deals with MIF and return more than only FORMAT fields
			done 
		fi
		#echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> tmp.merged.indels.tsv
	done < bedtools.multiintersection.2vcfs.indels.seu.slk.bed
	IFS=${old_IFS}
}

SEU_INDELS=$1
SLK_INDELS=$2
DUPLICATES_FILE=$3
getMergedINDELs "${SEU_INDELS}" "${SLK_INDELS}" "${DUPLICATES_FILE}"
exit

