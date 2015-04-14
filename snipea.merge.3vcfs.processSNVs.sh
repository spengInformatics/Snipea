#!/bin/bash
#USAGE: Initiate Snipea downstream processing starting from SNV

dir_script=`dirname $0`
source ${dir_script}/snipea.merge.3vcfs.functions.snvs.sh
source ${dir_script}/snipea.functions.static.sh
source ${dir_script}/snipea.merge.3vcfs.functions.indels.sh

function getMergedSNVs(){
	#Order will always be: SEURAT, STRELKA, MUTECT
	local SEU=$1
	local SLK=$2
	local MTC=$3
	local DUPLICATES_FILE=$4
	local TEMP_SNV_OUTFILE="tmp.merged.snvs.tsv"
	local CUR_LOCUS_FILE="CUR_LOCUS.SNVS.bed"
	local DUP_ALREADY_DONE_FILE="dup_snv_already_processed_as_dup.txt" ; >${DUP_ALREADY_DONE_FILE}

	grep -vE "^#"  ${SEU} | awk 'BEGIN {FS=OFS="\t"} {if($1 !~ /^#/ ) { $2=$2-1"\t"$2 ; print } else {print } }' | bedtools sort -i - > tmp.${SEU}.snvs.bed & 
	grep -vE "^#"  ${SLK} | awk 'BEGIN {FS=OFS="\t"} {if($1 !~ /^#/ ) { $2=$2-1"\t"$2 ; print } else {print } }' | bedtools sort -i - > tmp.${SLK}.snvs.bed &
	grep -vE "^#"  ${MTC} | awk 'BEGIN {FS=OFS="\t"} {if($1 !~ /^#/ ) { $2=$2-1"\t"$2 ; print } else {print } }' | bedtools sort -i - > tmp.${MTC}.snvs.bed &
	wait
	BTI_SEU=tmp.${SEU}.snvs.bed ; BTI_SLK=tmp.${SLK}.snvs.bed ; BTI_MTC=tmp.${MTC}.snvs.bed;
	bedtools multiinter -i ${BTI_SEU} ${BTI_SLK} ${BTI_MTC} | cut -f5 | sort | uniq -c  > countsOfCommonSNVsLoci.perSamples.bedtools.txt &
	bedtools multiinter -i ${BTI_SEU} ${BTI_SLK} ${BTI_MTC} | awk 'BEGIN{ FS=OFS="\t" } { if ($3-$2 <2) {print $0 } else { for(c=0;c<$3-$2;c++){ printf "%s\t%s\t%s",$1,$2+c,$2+c+1  ; for(j=4;j<=NF;j++) { printf "\t%s",$j } ; printf "\n" } } } ' > bedtools.multiintersection.3vcfs.snvs.seu.stlk.mtct.bed &
	wait

	old_IFS=$IFS      # save the field separator
	IFS=$'\t'
	new_IFS=${IFS}

	while ${IFS} read LINE ; do 
		echo -e "$LINE" > ${CUR_LOCUS_FILE} ;
		MCHR=`echo -e "${LINE}" | cut -f1`
		MPOS=`echo -e "${LINE}" | cut -f3`
		MID="."
		PRESENT_IN="`echo -e "$LINE" | cut -f 5 | sed "s/,/\t/g"`" ; ## \t : because we changed the IFS
		LOCUS="`echo -e "$LINE" | cut -f1-3`" ;
		
		if [[ `bedtools intersect -a ${CUR_LOCUS_FILE} -b ${DUP_ALREADY_DONE_FILE} -f 1 2>/dev/null | wc -l  ` -eq 1  ]] ; then continue ; fi
		if [[ -s ${DUPLICATES_FILE} ]] ## test if file is not empty,Has data 1, empty 0
		then
			if [[ `bedtools intersect -a ${CUR_LOCUS_FILE}  -b ${DUPLICATES_FILE} -f 1 | /usr/bin/wc -l ` -eq 1  ]] ;
			then 
				echo "WARNING: SAME LOCATION (${MCHR}:${MPOS}) FOUND in more than ONE line in VCF" 1>&2 ;
				## we removed the current position from the duplicated list and set it as Processed.
				bedtools intersect  -a ${DUPLICATES_FILE} -b ${CUR_LOCUS_FILE} >> ${DUP_ALREADY_DONE_FILE} &
				bedtools intersect  -a ${DUPLICATES_FILE} -b ${CUR_LOCUS_FILE} -v > TT ; mv TT ${DUPLICATES_FILE} &
				wait ; continue ## by continuing, this means we skip any data and info regarding this location
			fi
		fi ## end if of testing if duplicates file is not empty
		L1=() ; L2=(); L3=() ; CALLERS_COUNT=0; FLAG_TOOL=""
		for I in ${PRESENT_IN}
		do
		## get the LINE of Interest for each tool and add the LINE to an ARRAY; 
		##WARNING: THe format of the line has to be VCF (not BED); Do not forget to remove extra column position if bed format; e.g. ${SEU} should be in vcf format
		case ${I} in
			1) L1=(`bedtools intersect -a ${SEU} -b ${CUR_LOCUS_FILE}`)  ; CALLERS_COUNT=$((${CALLERS_COUNT}+1)) ; FLAG_TOOL=$(addToolnameAsFlag "${FLAG_TOOL[@]}" "SEURAT")  ;;
			2) L2=(`bedtools intersect -a ${SLK} -b ${CUR_LOCUS_FILE}`)  ; CALLERS_COUNT=$((${CALLERS_COUNT}+1)) ; FLAG_TOOL=$(addToolnameAsFlag "${FLAG_TOOL[@]}" "STRELKA") ;;
			3) L3=(`bedtools intersect -a ${MTC} -b ${CUR_LOCUS_FILE}`)  ; CALLERS_COUNT=$((${CALLERS_COUNT}+1)) ; FLAG_TOOL=$(addToolnameAsFlag "${FLAG_TOOL[@]}" "MUTECT")  ;;
			*) exit ;;
		esac
		done
		for I in ${PRESENT_IN}
		do
		## CAPTURE THE PREVALENCE for the FORMAT fields; PREVALENCE means we show the data first from SEURAT if exist, then from STRELKA if EXIST, then from MUTECT; AT LEAST one should exist
		case ${I} in
			1) if [[ "${L1[@]}" != "" ]] ; then PREVALENCE="SEU" ; fi ; break ;;
			2) if [[ "${L2[@]}" != "" ]] ; then PREVALENCE="SLK" ; fi ; break ;;
			3) if [[ "${L3[@]}" != "" ]] ; then PREVALENCE="MTC" ; fi ; break ;;
			*) exit ;;
		esac
		done
		MREF=$(getREF "${L1[3]}" "${L2[3]}" "${L3[3]}") 
		MALT=$(getALT "${L1[4]}" "${L2[4]}" "${L3[4]}")
		MQUAL=$(getQUAL "${L1[5]}" "${L2[7]}" "${L3[9]}" "${L3[10]}")
		MFILT="."
		MINFO=$(concatInfoFields "${L1[7]}" "${L2[7]}" "${L3[7]}")
		MFORMAT=$(getFORMAT "${L1[7]}" "${L2[7]}" "${L2[8]}" "${L2[9]}" "${L2[10]}" "${L3[7]}" "${L3[8]}" "${L3[10]}" "${L3[9]}" "${MREF}" "${MALT}" "${PREVALENCE}")


			if [[ ! ${MREF} =~ "," ]] ; 
			then 
				echo -e "${MCHR}\t${MPOS}\t${MID}\t${MREF}\t${MALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_SNV_OUTFILE}
			else
				echo -e "WARNING: this POSITION is AMBIGUOUS: ${MCHR}:${MPOS}\nFound called differently from one caller to another" 1>&2 ; 
				# splitting the POSITION into two lines: ## normally should not happen 
				local ARR_REF=(`echo ${MREF} | sed 's/,/\t/g'`) ; CREF=${#ARR_REF[@]} ; #echo -e "ARR_REF == ${ARR_REF[@]} \t CREF= ${CREF} ; \t ${MREF}" 1>&2 ;
				for (( L=1; L<=${CREF}; L++ ))
				do 
					local newMREF="$(echo ${MREF} | cut -d',' -f ${L})"
					local newMALT=`echo -e ${MALT} | cut -d',' -f ${L} `
					#echo -e "${MCHR}\t${MPOS}\t${MID}\t${newMREF}\t${newMALT}\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}"  1>&2
					echo -e "${MCHR}\t${MPOS}\t${MID}\t`echo -e \"${MREF}\" | cut -d',' -f${L}`\t`echo -e \"${MALT}\" | cut -d',' -f${L}`\t${MQUAL}\t${MFILT}\tSOMATIC;CALLERS_COUNT=${CALLERS_COUNT};${FLAG_TOOL};${MINFO};${MFORMAT}" >> ${TEMP_SNV_OUTFILE}
				done 
			fi
	done < bedtools.multiintersection.3vcfs.snvs.seu.stlk.mtct.bed
	IFS=${old_IFS}
}

SEU=$1
SLK=$2
MTC=$3
DUPLICATES_FILE=$4
getMergedSNVs "${SEU}" "${SLK}" "${MTC}" "${DUPLICATES_FILE}"
exit
