#!/bin/bash

SEU=$1
SLK=$2
MTC=$3
MYCURDIR=$4
PREFIX=$5
dir_scripts=$6
cd ${MYCURDIR}
FlagERR="OK"
for F in ${SEU} ${SLK} ${MTC}
do
	if [[ ! -e ${F} ]] ; then echo -e "FILE NOT FOUND: ${F}; ${PREFIX} snvs Venn Aborted." ; FlagERR="FNF" ; fi
done
if [[ ${FlagERR} == "FNF" ]] ; then exit -1 ; fi

cut -f1-2,4-5 ${SEU} | grep -vE "^#|^@" | sed 's/\t/_/g' > seu.venn
cut -f1-2,4-5 ${SLK} | grep -vE "^#|^@" | sed 's/\t/_/g' > slk.venn
cut -f1-2,4-5 ${MTC} | grep -vE "^#|^@" | sed 's/\t/_/g' > mtc.venn

INI_FILE=${MYCURDIR}/R_VennDiagram.ini
cat ${dir_scripts}/R_VennDiagram.ini  | sed -e 's,dir_infiles="",dir_infiles="'"${MYCURDIR}"'",'  > ${INI_FILE}
sed -i 's,prefix_output="",prefix_output="'"${PREFIX}"'.indels",' ${INI_FILE}
sed -i 's,comments2add="",comments2add="'"${PREFIX}"'.indels",' ${INI_FILE}

if [[ ${HOSTNAME} =~ "pnap" ]]
then
	/packages/R/2.15.2/bin/Rscript ${dir_scripts}/Venn.R ${INI_FILE}
fi
if [[ ${HOSTNAME} =~ "car" ]]
then
	/packages/R/2.15.2/bin/Rscript ${dir_scripts}/Venn.R ${INI_FILE}
fi

if [[ ${HOSTNAME} =~ "merckx" ]]
then
	## version R 3.0
	Rscript ${dir_scripts}/Venn.R ${INI_FILE}
else
	echo -e "ERROR : HOSTNAME NOT FOUND; Please ask you admin to add the host name to the script"
fi

