#!/bin/bash

function addHeader(){
	## input: $1 is SEURAT SNV file; $2 is STRELKA SNV file ; and $3 is MUTECT SNV
	## input $4 is the STRELKA INDEL; ONLY STRELKA because the header from seurat is already integrated in SNVS header; 
	## $5 is the command we used to run the merger (\$0 \$@ from the main script) 
	echo -e "$(headerSNVs $1 $2 $3 "$5")\n$(headerINDELs $4)\n$(addHeaderLine)"
}

function headerINDELs(){
	local SLK_INDELS_FILE=$1
	local HEADER_INDELS=`grep -E "^#" ${SLK_INDELS_FILE} | grep -vE "^##contig|^##fileformat|^##fileDate|^##source|^##startTime|^##reference|^#CHROM" | sed 's/<ID=/<ID=STRELKA_/g'`
	echo -e "${HEADER_INDELS}"
}

function headerSNVs(){
#input: $1 is SEURAT SNV file; $2 is STRELKA SNV file ; and $3 is MUTECT SNV
local SEU=$1 ; local SLK=$2 ; local MTC=$3 ; local command_run="$4"
local ALL=`cat ${SEU} ${SLK} ${MTC} | grep -E "^##|^#" | grep -vE "^##fileformat|^##fileDate|^#CHROM" | sort -u `
local contigs=`echo -e "${ALL}" | grep -E "^##contig"` # | sed 's/##contig/\n##contig/g'`
local filters=`echo -e "${ALL}" | grep -E "^##FILTER=" | sort -u` #  | sed 's/##FILTER/\n##FILTER/g'`
local ISEU=`cat ${SEU} | grep -E "^##" | grep -vE "^##contig|^##fileformat|^##fileDate|^#CHROM|^##FILTER" | sed 's/=<ID=/=<ID=SEURAT_/'`
local ISLK=`cat ${SLK} | grep -E "^##" | grep -vE "^##contig|^##fileformat|^##fileDate|^#CHROM|^##FILTER" | sed 's/=<ID=/=<ID=STRELKA_/'`
local IMTC=`cat ${MTC} | grep -E "^##" | grep -vE "^##contig|^##fileformat|^##fileDate|^#CHROM|^##FILTER" | sed 's/=<ID=/=<ID=MUTECT_/'`
ALL=`echo -e "${ISEU}\n${ISLK}\n${IMTC}"`
local formats=`echo -e "${ALL}" | grep -E "^##FORMAT" | sort ` #  | sed 's/##FORMAT/\n##FORMAT/g'`
local infos=`echo -e "${ALL}" | grep -E "##INFO" | sort` #  | sed 's/##INFO/\n##INFO/g'`
infos=$(correctNumberInSeuratField "${infos}")
local additionalInfos=`echo -e "${ALL}" | grep -vE "^##contig|^##FORMAT|^##INFO|^#CHROM" | sort -u` #  | sed 's/##/\n##/g'`


currentDate=`date +%Y%m%d`
echo -e "##fileformat=VCFv4.1
##fileDate=${currentDate}
##cmdLine_merger3vcfs=\"${command_run}\"
${additionalInfos}
${contigs}
${filters}
${infos}
$(customAddedInfoFlags)
${formats}
$(customAddedFormatFlags)"
}

function addHeaderLine(){
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR"
}

function customAddedInfoFlags(){
echo -e "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">
##INFO=<ID=CALLERS_COUNT,Number=1,Type=Integer,Description=\"Number of tools calling this variant event\">
##INFO=<ID=SEURAT,Number=0,Type=Flag,Description=\"Somatic mutation called by the SEURAT tool\">
##INFO=<ID=STRELKA,Number=0,Type=Flag,Description=\"Somatic mutation called by the STRELKA tool\">
##INFO=<ID=MUTECT,Number=0,Type=Flag,Description=\"Somatic mutation called by the MUTECT tool\">
##INFO=<ID=SEURAT_DP_NORMAL,Number=1,Type=Integer,Description=\"The depth of coverage in normal from SEURAT tool\">
##INFO=<ID=SEURAT_AR_NORMAL,Number=1,Type=Float,Description=\"Allele frequency of ALT allele in normal from SEURAT tool\">
##INFO=<ID=SEURAT_AD_NORMAL,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed in normal from SEURAT tool\">
##INFO=<ID=SEURAT_GT_NORMAL,Number=1,Type=String,Description=\"Genotype in normal from SEURAT tool\">
##INFO=<ID=SEURAT_DP_TUMOR,Number=1,Type=Integer,Description=\"The depth of coverage in tumor from SEURAT tool\">
##INFO=<ID=SEURAT_AR_TUMOR,Number=1,Type=Float,Description=\"Allele frequency of ALT allele in tumor from SEURAT tool\">
##INFO=<ID=SEURAT_AD_TUMOR,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed in tumor from SEURAT tool\">
##INFO=<ID=SEURAT_GT_TUMOR,Number=1,Type=String,Description=\"Genotype in tumor from SEURAT tool\">
##INFO=<ID=SEURAT_REF,Number=1,Type=String,Description=\"Reference INDEL captured by SEURAT tool\">
##INFO=<ID=SEURAT_ALT,Number=1,Type=String,Description=\"ALT INDEL captured by SEURAT tool\">
##INFO=<ID=STRELKA_DP_NORMAL,Number=1,Type=Integer,Description=\"The depth of coverage in normal from STRELKA tool\">
##INFO=<ID=STRELKA_AR_NORMAL,Number=1,Type=Float,Description=\"Allele frequency of ALT allele in normal from STRELKA tool\">
##INFO=<ID=STRELKA_AD_NORMAL,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed in normal from STRELKA tool\">
##INFO=<ID=STRELKA_GT_NORMAL,Number=1,Type=String,Description=\"Genotype in normal from STRELKA tool\">
##INFO=<ID=STRELKA_DP_TUMOR,Number=1,Type=Integer,Description=\"The depth of coverage in tumor from STRELKA tool\">
##INFO=<ID=STRELKA_AR_TUMOR,Number=1,Type=Float,Description=\"Allele frequency of ALT allele in tumor from STRELKA tool\">
##INFO=<ID=STRELKA_AD_TUMOR,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed in tumor from STRELKA tool\">
##INFO=<ID=STRELKA_GT_TUMOR,Number=1,Type=String,Description=\"Genotype in tumor from STRELKA tool\">
##INFO=<ID=MUTECT_DP_NORMAL,Number=1,Type=Integer,Description=\"The depth of coverage in normal from MUTECT tool\">
##INFO=<ID=MUTECT_AR_NORMAL,Number=1,Type=Float,Description=\"Allele frequency of ALT allele in normal from MUTECT tool\">
##INFO=<ID=MUTECT_AD_NORMAL,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed in normal from MUTECT tool\">
##INFO=<ID=MUTECT_GT_NORMAL,Number=1,Type=String,Description=\"Genotype in normal from MUTECT tool\">
##INFO=<ID=MUTECT_DP_TUMOR,Number=1,Type=Integer,Description=\"The depth of coverage in tumor from MUTECT tool\">
##INFO=<ID=MUTECT_AR_TUMOR,Number=1,Type=Float,Description=\"Allele frequency of ALT allele in tumor from MUTECT tool\">
##INFO=<ID=MUTECT_AD_TUMOR,Number=2,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed in tumor from MUTECT tool\">
##INFO=<ID=MUTECT_GT_TUMOR,Number=1,Type=String,Description=\"Genotype in tumor from MUTECT tool\">"

}

function customAddedFormatFlags(){
echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered) from chosen prevalent tool\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tier1 (used+filtered)\">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed from chosen prevalent tool\">
##FORMAT=<ID=AR,Number=1,Type=Float,Description=\"Allele frequency of ALT allele from chosen prevalent tool\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype from chosen prevalent tool\">"
}

function correctNumberInSeuratField(){
## input is lines with INFO field
local I="$1"
local ISEUnew=`echo -e "${I}" | sed -e 's/ALT_ALLELE_FORWARD_FRACTION,Number=4/ALT_ALLELE_FORWARD_FRACTION,Number=1/g ; s/ALT_ALLELE_FORWARD,Number=4/ALT_ALLELE_FORWARD,Number=1/g ; s/ALT_ALLELE_REVERSE_FRACTION,Number=4/ALT_ALLELE_REVERSE_FRACTION,Number=1/g ; s/ALT_ALLELE_REVERSE,Number=1/ALT_ALLELE_REVERSE,Number=1/g ; s/ALT_ALLELE_REVERSE,Number=4/ALT_ALLELE_REVERSE,Number=1/g ; s/ALT_ALLELE_TOTAL_FRACTION,Number=4/ALT_ALLELE_TOTAL_FRACTION,Number=1/g ; s/ALT_ALLELE_TOTAL,Number=4/ALT_ALLELE_TOTAL,Number=1/g ; s/REF_ALLELE_FORWARD,Number=4/REF_ALLELE_FORWARD,Number=1/g ; s/REF_ALLELE_REVERSE,Number=4/REF_ALLELE_REVERSE,Number=1/g ; s/REF_ALLELE_TOTAL,Number=4/REF_ALLELE_TOTAL,Number=1/g' `
echo -e "${ISEUnew}"

}

