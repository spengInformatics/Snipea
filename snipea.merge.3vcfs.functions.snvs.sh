#!/bin/bash

function addToolnameAsFlag(){
	## $1 id the LISTofFlags ; $2 is the TOOLNAME
	local FLAG_TOOL=$1; local TOOLNAME=$2
	if [[ ${FLAG_TOOL[@]} != "" ]] ; then FLAG_TOOL="${FLAG_TOOL[@]};${TOOLNAME}" ; else FLAG_TOOL="${TOOLNAME}" ; fi
	echo -e "${FLAG_TOOL[@]}"
}

## Add lines Captured to an ARRAY:
function addLINES(){
	local ARRAYLINES=(${1}) 
	local DATA="$2"
	if [[ ${ARRAYLINES[@]} == "" ]] ; then echo "YES" ; ARRAYLINES=( \"${DATA[@]}\" ) ; else ARRAYLINES=( ${ARRAYLINES[@]} \"${DATA[@]}\" ); fi
	if [[ ${#ARRAYLINES[@]} -gt 0 ]] ; then echo -e "${ARRAYLINES[@]}" ; else echo -e "" ; fi
}

function addAltString(){
	NSTRING="$1"
	DATA="$2"
#	echo -e "XXX ${NSTRING} XXX -lt DDD ${DATA} DDD ||||| ${#NSTRING} -lt ${#DATA} ;;;;;; `echo -e "${NSTRING}" | wc -c` -lt `echo -e "${DATA}" | wc -c`" 1>&2  

	if [[ $NSTRING == "" ]] ; then NSTRING="$DATA" ;
	elif [[ $NSTRING == $DATA ]] ; then NSTRING="$DATA"  ;
	else
		if [[ ${#DATA} -eq 1 ]]  ## here NSTRING (stands for NewString) greater than current ALT; we want to add, we need to check if the new ALT is already in the NSTRING
		then
			if [[ ! ${NSTRING[*]} =~ ${DATA} ]] ; then NSTRING="${NSTRING},${DATA}" ; fi
		elif [[ `echo ${NSTRING} | wc -c` -lt `echo ${DATA} | wc -c` ]]
		then
			old_IFS=$IFS ; IFS=$',' ; new_IFS=${IFS} ; # manage field separator
			local TMPALT="${NSTRING}" ; read -a ARRTMPALT <<< "${TMPALT}" ; 
			NSTRING=${DATA} ; DATA=${TMPALT}
			read -a ARRDATA <<< "${DATA}" ; read -a ARRNSTRING <<< "${NSTRING}";
			IFS=${old_IFS}
			for ITEM in ${DATA}
			do
				if [[ ! ${NSTRING[*]} =~ ${ITEM} ]] ; then NSTRING="${NSTRING},${ITEM}" ; fi
			done

		elif [[ `echo ${NSTRING} | wc -c` -gt `echo ${DATA} | wc -c` ]] 
		then
			old_IFS=$IFS ; IFS=$',' ; new_IFS=${IFS} ; # manage field separator
			read -a ARRDATA <<< "${DATA}" ; read -a ARRNSTRING <<< "${NSTRING}";
			IFS=${old_IFS}
			for ITEM in ${DATA}
			do
				if [[ ! ${NSTRING[*]} =~ ${ITEM} ]] ; then NSTRING="${NSTRING},${ITEM}" ; fi
			done
		elif [[ ${NSTRING} == ${DATA} ]] ; then NSTRING=${NSTRING};
		else
			old_IFS=$IFS ; IFS=$',' ; new_IFS=${IFS} ; # manage field separator
			read -a ARRDATA <<< "${DATA}" ; IFS=${old_IFS}
			for ITEM in ${DATA}
			do
				if [[ ! ${NSTRING[*]} =~ ${ITEM} ]] ; then NSTRING="${NSTRING},${ITEM}" ; fi
			done
		fi
	fi
	#NSTRING=`echo ${NSTRING[@]} | sed 's/\n/,/g'`
	echo -e "${NSTRING}"
}

function getREF(){
	local MREF="" # MALT stands for MergedREF fields
	if [[ $1 == $2 && $2 == $3 ]]
	then
		MREF=$1
	else 
		if [[ $1 != "" ]] 
		then
			IFSEU="`echo $1`"
			MREF=$(addAltString "${MREF[@]}" "${IFSEU[@]}")
		fi
		if [[ $2 != "" ]] 
		then
			IFSLK="`echo $2`"
			MREF=$(addAltString "${MREF[@]}" "${IFSLK[@]}")
		fi
		if [[ $3 != "" ]] 
		then
			IFMTC="`echo $3`"
			MREF=$(addAltString "${MREF[@]}" "${IFMTC[@]}")
		fi
	fi
	## WE WILL HAVE TO CHANGE THE CODE HERE TO PRINT ONLY THE UNIQUE ONE; WE WILL NEED TO LOOP OVER the array a to print its content
#	MREF=`echo ${MREF} | awk -F"," 'BEGIN {FS=OFS=","} {for(i=1;i<=NF;i=i+1){a[$i]=$i} ; if(length(a)>1){print $0} else {print $1}} ' | sed 's/ //g ; s/\n$//'` 
	echo -e "${MREF[@]}"
}

function getALT(){
	local MALT="" # MALT stands for MergedAltField
	if [[ $1 == $2 && $2 == $3 ]]
	then
		MALT=$1
	else 
		if [[ $1 != "" ]] 
		then
			IFSEU="`echo $1`"
			MALT=$(addAltString "${MALT[@]}" "${IFSEU[@]}")
		fi
		if [[ $2 != "" ]] 
		then
			IFSLK="`echo $2`"
			MALT=$(addAltString "${MALT[@]}" "${IFSLK[@]}")
		fi
		if [[ $3 != "" ]] 
		then
			IFMTC="`echo $3`"
			MALT=$(addAltString "${MALT[@]}" "${IFMTC[@]}")
		fi
	fi
	echo -e "${MALT[@]}"
}

function factorial() {
        printf '%d %s\n' $1 '1sA [q]sB [d1!<BdlA*sA1-lFx]sF lFx lAp' | dc
}

function getQUAL(){
	## Average quality score of all three callers for SNV (two for INDEL).
	local CALLERS_NUM=0
	local SEU_QUAL=0
	local SLK_QUAL=0
	local MTC_QUAL=0
	if [[ $1 != "" ]]  
	then 
		SEU_QUAL=$1 
		if [[ 1 -eq "$(echo "${SEU_QUAL} > 60" | bc -l)" ]]; then SEU_QUAL=60; fi
		CALLERS_NUM=$((${CALLERS_NUM}+1))
	fi
	if [[ $2 != "" ]]  
	then 
		local pos=`echo $2|awk -F";" '{for(i=1;i<=NF;i++){ if($i ~ /^QS._NT/) { print i ; break }}}'`
		SLK_QUAL=`echo $2|cut -d";" -f ${pos}|cut -d= -f2`
		if [[ 1 -eq "$(echo "${SLK_QUAL} > 60" | bc -l)" ]]; then SLK_QUAL=60; fi
		CALLERS_NUM=$((${CALLERS_NUM}+1))
	fi	
	if [[ $3 != "" ]]  
	then 
		local a=`echo $3|cut -d: -f2|cut -d, -f1`
		local c=`echo $3|cut -d: -f2|cut -d, -f2`
		local b=`echo $4|cut -d: -f2|cut -d, -f1`
		local d=`echo $4|cut -d: -f2|cut -d, -f2`
		local FISHER_TEST=`echo "$(factorial $(($a+$b)))/$(factorial $a)/$(factorial $b)*$(factorial $(($c+$d)))/$(factorial $c)/$(factorial $d)*$(factorial $(($a+$c)))*$(factorial $(($d+$b)))/$(factorial $(($a+$b+$c+$d)))"|bc -l`
		MTC_QUAL=`printf "%.2f" $(echo "-10*l(${FISHER_TEST})/l(10)"|bc -l)`
		if [[ 1 -eq "$(echo "${MTC_QUAL} > 60" | bc -l)" ]]; then MTC_QUAL=60; fi
		CALLERS_NUM=$((${CALLERS_NUM}+1))
	fi	
	local TOTAL=`echo "${SEU_QUAL}+${SLK_QUAL}+${MTC_QUAL}"|bc -l`
	QUAL=`printf "%.2f" $(echo "${TOTAL}/${CALLERS_NUM}"|bc -l)`
	echo -e "${QUAL}"
}

############
function getTAGvalueFromINFO(){
	## function to capture any TAG from VCF INFO field and return the Value of the TAG
	## $1 is the TAG (aka FLAG, aka FIELD, aka ID) to capture
	## $2 is the Tool VCF's INFO field
	local V=`echo $2 | awk -F";" -v TAG=$1 '{for(i=1;i<=NF;i++){ if($i ~ TAG) { print $i ; break }}}' | cut -d"=" -f2`
	echo -e "${V}"
}
function getTAGvalueFromFORMAT(){
	## function to capture any TAG from VCF FORMAT field and return the Value of the TAG
	## NOTE: this function is dedicated EXCLUSIVELY to VCF containing results from SOMATIC calls
	## THEREFORE: two columns should follow the FORMAT column: NORMAL and TUMOR columns
	## $1 is the TAG (aka FLAG, aka FIELD, aka ID) to capture
	## $2 is the Tool VCF's FORMAT column
	## $3 is the Sample column where the value has to be retrieved (here we do not care about the type of the sample)
	local TAG=$1 ; local FORMAT=$2 ; local SAMPLE=$3
	IDXTAG=`echo "${FORMAT}" | awk -F":" -v TAG=$TAG '{for(i=1;i<=NF;i++){ if($i ~ TAG) { print i ; break }}}'`
	TAGVALUE=`echo "${SAMPLE}" | cut -d":" -f${IDXTAG}`
	echo -e "${TAGVALUE}"
}

function getAD_fromSEU_fields_forTUMOR(){
	local TOTAL_REF=-1 ; local TOTAL_ALT=-1; local INFO=$1
	TOTREF=$( getTAGvalueFromINFO "DNA_REF_ALLELE_TOTAL" ${INFO} )
	TOTALT=$( getTAGvalueFromINFO "DNA_ALT_ALLELE_TOTAL" ${INFO} )
	echo -e "${TOTREF},${TOTALT}"
}

function getAD_fromSEU_fields_forNORMAL(){
	local TOTAL_REF=-1 ; local TOTAL_ALT=-1 ; local INFO=$1
	local DP1=$( getTAGvalueFromINFO "DP1" ${INFO} )
	local AR1=$( getTAGvalueFromINFO "AR1" ${INFO} )
	TOTALT=`echo ${DP1}*${AR1} | bc -l` ; TOTALT=`printf '%.0f\n' ${TOTALT}`
	TOTREF=$((${DP1}-${TOTALT}))
	echo -e "${TOTREF},${TOTALT}"
}

###########################################
function calculate_AR_forSLK(){
	## need the DP, AU, CU GU and TU fields from SLK  (this is in the FORMAT field all of them)
	## need the total number of ALT allele from SLK
	## need to know the REF from position (What About INDELs??? Here we only deal with SNV not indels from Strelka and Mutect)
	## so	## $1 is REF ; $2 is DP ; $3 is FORMAT ; $4 is SAMPLE_XX column
	local TOTAL_ALT="-1" ; local TOTAL_REF="-1"
	local REF="$1" ; local ALT="$2" ; local FORMAT="$3" ; local SAMPLE="$4" ;
	if [[ ${ALT} =~ "," ]] ; then ERRALT=${ALT} ; ALT=`echo ${ALT} | cut -d"," -f1` ; fi
	local TOTAL_REF=$(getNumberOfALLELESfromSLK "$REF" "$FORMAT" "$SAMPLE")
	local TOTAL_ALT=$(getNumberOfALLELESfromSLK "$ALT" "$FORMAT" "$SAMPLE")
	if [[ ${TOTAL_ALT} == "-1" || ${TOTAL_REF} == "-1" ]] ; then echo -e "ERROR in calculate_AR_forSLK ; return(ed) value is -1\nALT was: ${ERRALT}" >> logs.snipea.merge.3vcfxs.functions.sh.log ; echo -e "-1"; exit ; fi
	RESULT="`echo -e "scale=3; ${TOTAL_ALT}/(${TOTAL_ALT}+${TOTAL_REF})" | bc `"
	AR=`printf '%.3f\n' ${RESULT}` ; if [[ $? -ne 0 ]] ; then echo -e "-1" ; exit ; fi
	echo -e "${AR}"
}

function getNumberOfALLELESfromSLK(){
	local TOTAL_ALT=-1
	local ALLELE=$1 ; local FORMAT="$2" ; local SAMPLE="$3"
	## AU, CU, GU and TU contain two numbers separated with a comma (tier1,tier2); We have chosen to keep tier1; the most stringent tier.
	if [[ ${ALT} =~ "," ]] ; then ERRALT=${ALT} ; ALT=`echo ${ALT} | cut -d"," -f1` ; fi
	case $ALLELE in
		A) AU=$(getTAGvalueFromFORMAT "AU" "${FORMAT}" "${SAMPLE}") ; AU=`echo ${AU} | cut -d"," -f1 ` ; TOTALT=${AU} ;;
		C) CU=$(getTAGvalueFromFORMAT "CU" "${FORMAT}" "${SAMPLE}") ; CU=`echo ${CU} | cut -d"," -f1 ` ; TOTALT=${CU} ;;
		G) GU=$(getTAGvalueFromFORMAT "GU" "${FORMAT}" "${SAMPLE}") ; GU=`echo ${GU} | cut -d"," -f1 ` ; TOTALT=${GU} ;;
		T) TU=$(getTAGvalueFromFORMAT "TU" "${FORMAT}" "${SAMPLE}") ; TU=`echo ${TU} | cut -d"," -f1 ` ; TOTALT=${TU} ;;
		*) TOTALT="-1" ;;
	esac
	echo -e "${TOTALT}"
}

function calculate_GT_using_AR_value(){
	## need AR for appropriate sample; in Seurat (this function is SEURAT's SPECIFIC, everything is in INFO field
	local AR=$1 ; ##TODO : Check if both are numbers
	local thld_min=0.200; local thld_max=0.900;
	if [[ `echo ${AR} '>' ${thld_max} | bc -l ` -eq 1 || `echo ${AR} '==' ${thld_max} | bc -l ` -eq 1 ]] ; then echo -e "1/1" ; 
	elif [[ `echo ${AR} '<' ${thld_max} | bc -l ` -eq 1 && `echo ${AR} '>' ${thld_min} | bc -l ` -eq 1 ]] ; then echo -e "0/1" ; 
	elif [[ `echo ${AR} '<' ${thld_min} | bc -l ` -eq 1 || `echo ${AR} '==' ${thld_min} | bc -l ` -eq 1 ]] ; then echo -e "0/0" ; 
	else echo -e "ERROR_WITH_AR_VALUE___GT_NOT_CAPTURED" ;
	fi
}

## functions for the snipea.merge.3vcfs.sh script files
function addToMIF(){
	MIF="$1"
	shift
	for DATA in $* ; do if [[ $MIF == "" ]] ; then MIF="$DATA" ; else MIF="${MIF[@]};${DATA[@]}" ; fi; done
	echo ${MIF[@]}
}

function concatInfoFields(){
	#$1 is INFO field from Seurat, #$2 is INFO field from Strelka, #$1 is INFO field from Mutect, 
	local MIF="" # MIF stdands for MergedInfoField
	if [[ $1 != "" ]] 
	then
		IFSEU="`echo $1 | sed 's/;$// ; s/^/SEURAT_/ ; s/;/;SEURAT_/g'`"
		MIF=$(addToMIF "${MIF[@]}" "${IFSEU[@]}")
	fi
	if [[ $2 != "" ]] 
	then
		IFSLK="`echo $2 | sed 's/;$// ; s/^/STRELKA_/ ; s/;/;STRELKA_/g'`"
		MIF=$(addToMIF "${MIF[@]}" "${IFSLK[@]}")
	fi
	if [[ $3 != "" ]] 
	then
		IFMTC="`echo $3 | sed 's/;$// ; s/^/MUTECT_/ ; s/;/;MUTECT_/g'`"
		MIF=$(addToMIF "${MIF[@]}" "${IFMTC[@]}")
	fi
#	MIF="`echo ${MIF[@]} | sed 's/;$//'`"
	echo -e "${MIF[@]}"
}


function getFORMAT() {
	## $1 is INFO field from Seurat
	## $2 is INFO field from Strelka
	## $3 is FORMAT field from Strelka
	## $4 is NORMAL field from Strelka
	## $5 is TUMOR field from Strelka
	## $6 is INFO field from Mutect
	## $7 is FORMAT field from Mutect
	## $8 is NORMAL field from Mutect
	## $9 is TUMOR field from Mutect
	## $10 is the REF 
	## $11 is the ALT
	## $12 is the PREVALENCE
	#echo -e "${1}\t${2}\t${3}\t${4}\t${5}\t${6}\t${7}\t${8}\t${9}\t${10}\t${11}\t${12}"
	local REF=${10}
	local ALT=${11}
	local PREVALENCE=${12}
	local MIF=""

	if [[ $1 != "" ]]
	then
			DPSEUnormal=$(getTAGvalueFromINFO "DP1" ${1} ) ; DPSEUtumor=$(getTAGvalueFromINFO "DP2" ${1} ) ; 
			AR1SEUval=$(getTAGvalueFromINFO "AR1" ${1}) ; AR2SEUval=$(getTAGvalueFromINFO "AR2" ${1}) ; 
			ADSEUnormal="`echo "scale=0; ${AR1SEUval}*${DPSEUnormal}" | bc `"; ADSEUnormal=`printf '%.0f\n' ${ADSEUnormal}` ; 
			ADSEUnormal="$( getAD_fromSEU_fields_forNORMAL ${1} )" 
			ADSEUtumor="$( getAD_fromSEU_fields_forTUMOR ${1} )"
			GTSEUnormal=$(calculate_GT_using_AR_value ${AR1SEUval}) ; GTSEUtumor=$(calculate_GT_using_AR_value ${AR2SEUval}) ; 
			MIF=$(addToMIF "${MIF[@]}" "SEURAT_DP_NORMAL=${DPSEUnormal}" "SEURAT_DP_TUMOR=${DPSEUtumor}" "SEURAT_AR_NORMAL=${AR1SEUval}" "SEURAT_AR_TUMOR=${AR2SEUval}" "SEURAT_AD_NORMAL=${ADSEUnormal}" "SEURAT_AD_TUMOR=${ADSEUtumor}" "SEURAT_GT_NORMAL=${GTSEUnormal}" "SEURAT_GT_TUMOR=${GTSEUtumor}" )
	fi
	if [[ $2 != "" ]] 
	then

			DPSLKnormal=$(getTAGvalueFromFORMAT "DP" ${3} ${4}) ; DPSLKtumor=$(getTAGvalueFromFORMAT "DP" ${3} ${5}) ; local DP_NORMAL=${DPSLKnormal}; 
			ARSLKnormal=$(calculate_AR_forSLK "${REF}" "${ALT}" "${3}" "${4}") ; ARSLKtumor=$(calculate_AR_forSLK "${REF}" "${ALT}" "${3}" "${5}") ; 
			ADSLKnormal="$(getNumberOfALLELESfromSLK ${REF} ${3} ${4} ),$(getNumberOfALLELESfromSLK  ${ALT} ${3} ${4} )"
			ADSLKtumor="$(getNumberOfALLELESfromSLK ${REF} ${3} ${5} ),$(getNumberOfALLELESfromSLK  ${ALT} ${3} ${5} )"
			GTSLKnormal=$(calculate_GT_using_AR_value ${ARSLKnormal}) ; GTSLKtumor=$(calculate_GT_using_AR_value ${ARSLKtumor}) ; 
			MIF=$(addToMIF "${MIF[@]}" "STRELKA_DP_NORMAL=${DPSLKnormal}" "STRELKA_DP_TUMOR=${DPSLKtumor}" "STRELKA_AR_NORMAL=${ARSLKnormal}" "STRELKA_AR_TUMOR=${ARSLKtumor}" "STRELKA_AD_NORMAL=${ADSLKnormal}" "STRELKA_AD_TUMOR=${ADSLKtumor}" "STRELKA_GT_NORMAL=${GTSLKnormal}" "STRELKA_GT_TUMOR=${GTSLKtumor}"  )
	fi
	if [[ $6 != "" ]]
	then
			DPMTCnormal=$(getTAGvalueFromFORMAT "DP" ${7} ${8}) ; DPMTCtumor=$(getTAGvalueFromFORMAT "DP" ${7} ${9}) ; 
			ARMTCnormal=$(getTAGvalueFromFORMAT "FA" ${7} ${8}) ; ARMTCtumor=$(getTAGvalueFromFORMAT "FA" ${7} ${9}) ; 
			ADMTCnormal=`echo $(getTAGvalueFromFORMAT "AD" ${7} ${8})` ; ADMTCtumor=`echo $(getTAGvalueFromFORMAT "AD" ${7} ${9})` ; local AD_NORMAL=${ADMTCnormal} ; local AD_TUMOR=${ADMTCtumor}
			GTMTCnormal=$(getTAGvalueFromFORMAT "GT" ${7} ${8}) ; GTMTCtumor=$(getTAGvalueFromFORMAT "GT" ${7} ${9}) ;  
			MIF=$(addToMIF "${MIF[@]}" "MUTECT_DP_NORMAL=${DPMTCnormal}" "MUTECT_DP_TUMOR=${DPMTCtumor}" "MUTECT_AR_NORMAL=${ARMTCnormal}" "MUTECT_AR_TUMOR=${ARMTCtumor}" "MUTECT_AD_NORMAL=${ADMTCnormal}" "MUTECT_AD_TUMOR=${ADMTCtumor}" "MUTECT_GT_NORMAL=${GTMTCnormal}" "MUTECT_GT_TUMOR=${GTMTCtumor}" )
	fi

	case ${PREVALENCE} in
		SEU) 
			local DP_NORMAL=${DPSEUnormal}; local DP_TUMOR="${DPSEUtumor}"
			local AR_NORMAL=${AR1SEUval} ; local AR_TUMOR=${AR2SEUval}
			local AD_NORMAL=${ADSEUnormal} ; local AD_TUMOR=${ADSEUtumor} ; 
			local GT_NORMAL="${GTSEUnormal}" ;  local GT_TUMOR="${GTSEUtumor}" ;
			;;
		SLK) 
			local DP_NORMAL=${DPSLKnormal}; local DP_TUMOR="${DPSLKtumor}"
			local AR_NORMAL=${ARSLKnormal} ; local AR_TUMOR=${ARSLKtumor}
			local AD_NORMAL=${ADSLKnormal} ; local AD_TUMOR=${ADSLKtumor} ; 
			local GT_NORMAL="${GTSLKnormal}" ; local GT_TUMOR="${GTSLKtumor}" ;
			;;
		MTC) 
			local DP_NORMAL=${DPMTCnormal}; local DP_TUMOR=${DPMTCtumor}
			local AR_NORMAL=${ARMTCnormal} ; local AR_TUMOR=${ARMTCtumor}
			local AD_NORMAL=${ADMTCnormal} ; local AD_TUMOR=${ADMTCtumor}
			local GT_NORMAL="${GTMTCnormal}" ; local GT_TUMOR="${GTMTCtumor}" ; 
			;;
		*) echo -e "ERROR_CAPTURING_FORMAT_FIELDS; Check with your System Administrator or Pipeline Administrator" ; exit ;;
	esac
	
	echo -e "${MIF}\tGT:DP:AR:AD\t${GT_NORMAL}:${DP_NORMAL}:${AR_NORMAL}:${AD_NORMAL}\t${GT_TUMOR}:${DP_TUMOR}:${AR_TUMOR}:${AD_TUMOR}"
	exit
}

