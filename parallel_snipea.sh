#!/bin/bash

Location_DIR=/scratch/speng/projects/SNV_callers/trial/Miseq

SCRIPTS='/scratch/speng/projects/SNV_callers/trial/jobScript'
COUNT=0

for sample in `find ${Location_DIR}  -name "*.seurat.vcf"`
do
       let COUNT=$COUNT+1
	   dirName=`dirname $sample`
	   fileName=`basename $sample|cut -d. -f1`
       Seurat_Basename=$fileName
#	   cd $dirName
	     echo $dirName
     qsub -v SEURAT_BASENAME=$Seurat_Basename,DIR_NAME=$dirName $SCRIPTS/Merge_VCFs.pbs
done

echo $COUNT
