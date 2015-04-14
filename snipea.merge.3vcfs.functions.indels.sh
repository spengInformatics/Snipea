#!/bin/bash

function indels_getREF(){
	local MREF="" # MREF stands for MergedREF field
	if [[ $1 == $2 ]]
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
	## WE WILL HAVE TO CHANGE THE CODE HERE TO PRINT ONLY THE UNIQUE ONE; WE WILL NEED TO LOOP OVER the array a to print its content; Normally we do not need this following line anymore; Test have to be done to see if we can get rid of it. ##TODO
	MREF=`echo ${MREF} | awk -F"," 'BEGIN {FS=OFS=","} {for(i=1;i<=NF;i=i+1){a[$i]=$i} ; if(length(a)>1){print $0} else {print $1}} ' | sed 's/ //g'` 
	echo -e "${MREF[@]}"
}

function indels_getALT(){
	local MALT="" # MALT stands for MergedAlt field
	if [[ $1 == $2 ]]
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
	## remove redundant information from if all ALT are same value; Otherwise keeps all.
	MALT=`echo ${MALT} | awk -F"," 'BEGIN {FS=OFS=","} {for(i=1;i<=NF;i=i+1){a[$i]=$i} ; if(length(a)>1){print $0}else {print $1}} '`
	echo -e "${MALT[@]}"
}


function indels_getTAGvalueFromFORMAT(){
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

function indels_calculate_AR_forSLK(){
	## need the TAR and TIR fields from SLK  (this is in the FORMAT field all of them)
	## need the total number of ALT allele from SLK
	## so	## $1 is FORMAT ; $2 is SAMPLE_XX column
	local FORMAT="$1" ; local SAMPLE="$2" ;
	local TOTAL_ALT="-1" ; local TOTAL_REF="-1"
	local TOTAL_REF=`echo -e "$(indels_getTAGvalueFromFORMAT "TAR" "$FORMAT" "$SAMPLE")" | cut -d"," -f1`
	local TOTAL_ALT=`echo -e "$(indels_getTAGvalueFromFORMAT "TIR" "$FORMAT" "$SAMPLE")" | cut -d"," -f1`
	if [[ ${ALT} =~ "," ]] ; then ERRALT=${ALT} ; ALT=`echo ${ALT} | cut -d"," -f1` ; fi
	if [[ ${TOTAL_ALT} == "-1" || ${TOTAL_REF} == "-1" ]] ; then echo -e "ERROR in calculate_AR_forSLK ; return(ed) value is -1 " >> logs.snipea.merge.3vcfs.functions.sh.log ; echo -e "-1"; exit ; fi
	if [[ `echo ${TOTAL_REF} '==' 0 | bc -l ` -eq 1  &&  `echo ${TOTAL_ALT} '==' 0 | bc -l ` -eq 1 ]] ; then echo -e "0.000" ; exit ; fi ## we need this line because in the NORMAL sample, it may happen that both TAR and TIR have the value of ZERO
	RESULT="`echo -e "scale=3; ${TOTAL_ALT}/(${TOTAL_ALT}+${TOTAL_REF})" | bc `"
	AR=`printf '%.3f\n' ${RESULT}` ; if [[ $? -ne 0 ]] ; then echo -e "-1" ; exit ; fi
	echo -e "${AR}"

}


function indels_calculate_GT_using_AR_value(){
	## need AR for appropriate sample; in Seurat (this functio nis SEURAT's SPECIFIC, everything is in INFO field
	local AR=$1 ; ##TODO : Check if both are numbers
	local thld_min=0.200; local thld_max=0.900;
	if [[ `echo ${AR} '>' ${thld_max} | bc -l ` -eq 1 || `echo ${AR} '==' ${thld_max} | bc -l ` -eq 1 ]] ; then echo -e "1/1" ; 
	elif [[ `echo ${AR} '<' ${thld_max} | bc -l ` -eq 1 && `echo ${AR} '>' ${thld_min} | bc -l ` -eq 1 ]] ; then echo -e "0/1" ; 
	elif [[ `echo ${AR} '<' ${thld_min} | bc -l ` -eq 1 || `echo ${AR} '==' ${thld_min} | bc -l ` -eq 1 ]] ; then echo -e "0/0" ; 
	else echo -e "./." ;
	fi
}


function indels_getFORMAT() {
	## $1 is INFO field from Seurat
	## $2 is INFO field from Strelka
	## $3 is FORMAT field from Strelka
	## $4 is NORMAL field from Strelka
	## $5 is TUMOR field from Strelka
	## $6 is the REF 
	## $7 is the ALT
	## $8 is the PREVALENCE
	#echo -e "${1}\t${2}\t${3}\t${4}\t${5}\t${6}\t${7}\t${8}\t${9}\t${10}\t${11}\t${12}"
	local REF=${6}
	local ALT=${7}
	local PREVALENCE=${8}
	local MIF=""
	local FLAG_INDEL_SEU=FALSE

	if [[ "$1" != "" ]]
	then
			if [[ "$1" =~ "AR1=" ]]
			then
				DPSEUnormal=$(getTAGvalueFromINFO "DP1" ${1} ) ; DPSEUtumor=$(getTAGvalueFromINFO "DP2" ${1} ) ; 
				AR1SEUval=$(getTAGvalueFromINFO "AR1" ${1}) ; AR2SEUval=$(getTAGvalueFromINFO "AR2" ${1}) ; 
				ADSEUnormal="$( getAD_fromSEU_fields_forNORMAL ${1} )" 
				ADSEUtumor="$( getAD_fromSEU_fields_forTUMOR ${1} )"
				GTSEUnormal=$(indels_calculate_GT_using_AR_value ${AR1SEUval}) ; GTSEUtumor=$(indels_calculate_GT_using_AR_value ${AR2SEUval}) ; 
				MIF=$(addToMIF "${MIF[@]}" "SEURAT_DP_NORMAL=${DPSEUnormal}" "SEURAT_DP_TUMOR=${DPSEUtumor}" "SEURAT_AR_NORMAL=${AR1SEUval}" "SEURAT_AR_TUMOR=${AR2SEUval}" "SEURAT_AD_NORMAL=${ADSEUnormal}" "SEURAT_AD_TUMOR=${ADSEUtumor}" "SEURAT_GT_NORMAL=${GTSEUnormal}" "SEURAT_GT_TUMOR=${GTSEUtumor}" )
			else
				FLAG_INDEL_SEU=TRUE
			fi
			if [[ ${FLAG_INDEL_SEU} == "TRUE" &&  $2 == "" ]] ; then echo -e "\tGT:DP:AR:AD\t./.:.:.:.\t./.:.:.:.\n" ; exit
			elif [[ ${PREVALENCE} == "SEURAT_DUPLICATE" ]] ; then PREVALENCE="SEU" ;
			elif [[ ${FLAG_INDEL_SEU}=="TRUE" &&  $2 != "" ]] ; then PREVALENCE="SLK" ;
			fi
	fi

	if [[ $2 != "" ]] 
	then
			local TAR_NORMAL=`echo $(indels_getTAGvalueFromFORMAT "TAR" ${3} ${4}) | cut -d"," -f1`
			local TAR_TUMOR=`echo $(indels_getTAGvalueFromFORMAT "TAR" ${3} ${5}) | cut -d"," -f1`
			local TIR_NORMAL=`echo $(indels_getTAGvalueFromFORMAT "TIR" ${3} ${4}) | cut -d"," -f1`
			local TIR_TUMOR=`echo $(indels_getTAGvalueFromFORMAT "TIR" ${3} ${5}) | cut -d"," -f1`
			DPSLKnormal=$(indels_getTAGvalueFromFORMAT "DP" ${3} ${4}) ; DPSLKtumor=$(indels_getTAGvalueFromFORMAT "DP" ${3} ${5})
			if [[ ${TAR_NORMAL} -gt ${DPSLKnormal} ]] ; then DPSLKnormal=${TAR_NORMAL} ; fi
			if [[ ${TAR_TUMOR} -gt ${DPSLKtumor} ]] ; then DPSLKtumor=${TAR_TUMOR} ; fi
			ARSLKnormal=$(indels_calculate_AR_forSLK "${3}" "${4}") ; ARSLKtumor=$(indels_calculate_AR_forSLK "${3}" "${5}") ; 
			ADSLKnormal="${TAR_NORMAL},${TIR_NORMAL}" ; ADSLKtumor="${TAR_TUMOR},${TIR_TUMOR}"
			GTSLKnormal=$(indels_calculate_GT_using_AR_value ${ARSLKnormal}) ; GTSLKtumor=$(indels_calculate_GT_using_AR_value ${ARSLKtumor}) ; 
			MIF=$(addToMIF "${MIF[@]}" "STRELKA_DP_NORMAL=${DPSLKnormal}" "STRELKA_DP_TUMOR=${DPSLKtumor}" "STRELKA_AR_NORMAL=${ARSLKnormal}" "STRELKA_AR_TUMOR=${ARSLKtumor}" "STRELKA_AD_NORMAL=${ADSLKnormal}" "STRELKA_AD_TUMOR=${ADSLKtumor}" "STRELKA_GT_NORMAL=${GTSLKnormal}" "STRELKA_GT_TUMOR=${GTSLKtumor}"  )
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
		*) echo -e "ERROR_CAPTURING_FORMAT_FIELDS; Check with your System Administrator or Pipeline Administrator" ; exit ;;
	esac
		if [[ ${MIF[@]} != "" ]] ; then MIF=";${MIF}" ; fi ; ## allows to keep the delimiter between MIF features
		echo -e "${MIF}\tGT:DP:AR:AD\t${GT_NORMAL}:${DP_NORMAL}:${AR_NORMAL}:${AD_NORMAL}\t${GT_TUMOR}:${DP_TUMOR}:${AR_TUMOR}:${AD_TUMOR}\n"

}


