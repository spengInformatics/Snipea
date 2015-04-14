#!/bin/bash
dir_script="";
keepVennInputFiles="no"

function usage(){
	echo "USAGE:"
	echo -e "/Snipea.main.sh --seusnv Seurat.filtered.SNVs.vcf  --slksnv Strelka.filtered.SNVs.vcf --mtcsnv Mutect.SNVs.vcf  --seuindel  Seurat.filtered.INDELs.vcf --slkindel Strelka.filtered.INDELs.vcf --dirscript \"${0}\""
	echo -e "OR (both versions of the run are mutually exclusive. This version uses a seurat's file containing both snvs and indels within the file"
	echo -e "/Snipea.main.sh --seu Seurat.filtered.snvs_and_indelsvcf  --slksnv Strelka.filtered.SNVs.vcf --mtcsnv Mutect.SNVs.vcf  --slkindel Strelka.filtered.INDELs.vcf --dirscript \"${0}\""
	echo -e "\nPRE-REQUISITE: We assume that the following MANDATORY files are provided:"
	echo -e "Seurat_SNVs\nStrelka_SNVs\nMutect_SNVs\nSeurat_INDELs\nStrelka_INDELs"
	echo -e "WARNING: If one of these files is not provided the script will FAIL ... or will output with wrong results/data."
	echo -e "Other options available:\n--keepVennInputFiles : possible value: yes or no ; if yes, the input file for the Venn R script will not be deleted\n"
}

function getOptions(){
	# options may be followed by one colon to indicate they have a required argument
	if ! options=`getopt -o h: -l seu:,seusnv:,slksnv:,mtcsnv:,seuindel:,slkindel:,outprefix:,dirscript:,keepVennInputFiles: -- "$@"`
	then
	# something went wrong, getopt will put out an error message for us
		echo "ERROR in Arguments" ; usage
		exit -1
	fi
	eval set -- "$options"
	while [[ $# -gt 0 ]]
	do
		# for options with required arguments, an additional shift is required
		case $1 in
		--seu) SEU=$2 ;  shift;;
		--seusnv) SEU_SNV=$2 ;  shift;;
		--slksnv) SLK_SNV=$2 ;  shift ;;
		--mtcsnv) MTC_SNV=$2 ;  shift ;;
		--seuindel) SEU_INDEL=$2 ;  shift ;;
		--slkindel) SLK_INDEL="$2" ;  shift ;;
		--outprefix) PREFIX="$2" ; shift ;;
		--dirscript) dir_script="`echo $2 | sed 's/\/$//'`" ; shift ;;
		--keepVennInputFiles) keepVennInputFiles="`echo $2 | sed 's/\/$//' | tr '[A-Z]' 'a-z'`" ; shift ;; ##no is the default 
		-h) usage ; exit ; shift ;;
		--help) usage ; exit ; shift ;;
		(--) shift ; echo "--" ;;
		(-*) echo "$0: error - unrecognized option $1" 1>&2; exit -1 ;;
		(*) break ; echo "$0: error --- unrecognized option $1" 1>&2 ;;
		esac
		shift
	done
	#input recap
	echo -e "INPUTS:\nseu\t\t${SEU}\nseusnv\t\t${SEU_SNV}\nslksnv\t\t${SLK_SNV}\nmtcsnv\t\t${MTC_SNV}\nseuindel\t${SEU_INDEL}\nslkindel\t${SLK_INDEL}\n"
}


###############
#### MAIN #####
###############
##capture the options for this script
getOptions $@
#########################
### PRE-PROCESSING ######
#########################
if [[ ${dir_script} == "" ]] ; then echo -e "\n*********************\nERROR: DIR SCRIPT MISSING; PLEASE provide  << dirscript>> option. Aborting\n*********************" ; exit -1 ; elif [[ ! -e ${dir_script} ]] ; then echo -e "\n*********************\nDIR SCRIPT NOT FOUND << ${dir_script}\n********************* >> " ; exit -1 ;fi
export dir_script=${dir_script}
source ${dir_script}/snipea.functions.static.sh 
source ${dir_script}/snipea.merge.3vcfs.functions.main.sh

command_run="$0 $@"
#remove any previous existed tmp files 
rm tmp.merged.snvs.tsv tmp.merged.indels.tsv duplicated_positions.snvs.bed duplicated_positions.indels.bed > /dev/null 2>&1
if [[ ${SEU} != "" && -e ${SEU} ]] ; then  splitSNVsINDELsFromSeuratFile ${SEU} ; fi
if [[ "$(checkFile ${SEU_SNV} ${SLK_SNV} ${MTC_SNV} ${SEU_INDEL} ${SLK_INDEL})" == "FNF" ]] ; then echo -e "INPUT FILE(S) NOT FOUND. Aborting."; exit -1; fi
##filtering Mutect ; Removing the REJECT calls
echo -e "Removing REJECT from MuTecT vcf file ..." 1>&2 
grep -vE "REJECT" ${MTC_SNV} > `basename ${MTC_SNV} "_MuTect_All.vcf"`.Mutect.filt.vcf
MTC_SNV=`basename ${MTC_SNV} "_MuTect_All.vcf"`.Mutect.filt.vcf
echo ${MTC_SNV}
echo -e "${listDupPosInSnvs}\n${listDupPosInSnvs}" > duplicated_positions.both_snvs_and_indels.bed 
captureDuplicatedPositionsSnvs ${SEU_SNV} ${SLK_SNV} ${MTC_SNV} 	## create the following file: duplicated_positions.snvs.bed 
captureDuplicatedPositionsIndels ${SEU_INDEL} ${SLK_INDEL}		## create the following file: duplicated_positions.indels.bed
cat duplicated_positions.snvs.bed duplicated_positions.indels.bed | sort > duplicated_positions.both_snvs_and_indels.bed 

PREFIX=$(getPrefixOutput "${PREFIX}" "${SEU_SNV}")
echo -e "counts of snvs/indels:"
countEventsPerFile ${PREFIX} ${SEU_SNV} ${SLK_SNV} ${MTC_SNV} ${SEU_INDEL} ${SLK_INDEL}

#####################
### PROCESSING  #####
#####################
echo "processing SNVs in background..."
bash ${dir_script}/snipea.merge.3vcfs.processSNVs.sh "${SEU_SNV}" "${SLK_SNV}" "${MTC_SNV}" "duplicated_positions.snvs.bed"  &
echo "processing INDELs in background ..."
bash ${dir_script}/snipea.merge.3vcfs.processINDELs.sh "${SEU_INDEL}" "${SLK_INDEL}" "duplicated_positions.indels.bed"  &
echo "creating Venn Diagrams in background ..."
bash ${dir_script}/snipea.runGetVenn.snvs.sh "${SEU_SNV}" "${SLK_SNV}" "${MTC_SNV}" "`pwd`" "${PREFIX}" "${dir_script}" >/dev/null 2>&1 &
bash ${dir_script}/snipea.runGetVenn.indels.sh "${SEU_INDEL}" "${SLK_INDEL}" "`pwd`" "${PREFIX}" "${dir_script}" >/dev/null 2>&1 &
echo -e ">>> wait ... work in progress ... "
wait
##########################
### POST-PROCESSING ######
##########################
echo -e "creating header ..."
echo -e "$(addHeader ${SEU_SNV} ${SLK_SNV} ${MTC_SNV} ${SLK_INDEL} "${command_run}")" > header.vcf
echo -e "concatenating data and saving merged vcf file ..."
cat header.vcf  tmp.merged.snvs.tsv  tmp.merged.indels.tsv  | sed '/^$/d' >  "${PREFIX}.merge.vcf"
echo -e "reorganizing files ..."
mkdir Venns >/dev/null 2>&1 ; mv *.png *.pdf *L.tsv -t Venns >/dev/null 2>&1
echo -e "cleaning up temp files ..."
#cleaningTempFiles ${keepVennInputFiles}
exit
