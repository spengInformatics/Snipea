#!/bin/bash
#####################################################################
##### Author: Sen Peng #####
##### Date: December 2014#####
#######################################
#@@@@@@@@@@@@@@@@#
## INTRODUCTION ##
#@@@@@@@@@@@@@@@@#

Snipea (SNv Integration, Prioritization, Ensemble, and Annotation) is a tool to further identify, prioritize and annotate somatic SNV from patient tumor (Other Biological samples). Our tool will apply an ensemble-based approach to improve the accuracy and annotation of variant callers. “Snipea” calculates the authenticity and statistical significance that an SNV is causative for a query disease and hence provides a means of prioritizing candidate SNVs. The ranking algorithm will take into account the callers_agreement, clinical_impact (non-synonymous SNV), quality score from each tool and public_database_status.

#@@@@@@@@@@@@@@@@#
## REQUIREMENTS ##
#@@@@@@@@@@@@@@@@#
This software requires the following tools/software to be installed

bedtools version 2.17 or up
linux commands: cat, case, cd , cut, sed, grep, awk, sort, uniq, getopt, wc, mv, cp, head , pwd, echo, for while, if, bash, read, mkdir, wait, do, rm, local.
R version 2.15.2 with features TIFF, PNG and CAIRO enabled; NOTE: R version 3.0 has not been tested but should work.
R packages for version 2.15.2 of R: 
		VennDiagram 1.6.5 or up
		Vennerable 2.2 or up
		gplots 2.11.0 or up
		optparse --> run  install.packages("optparse") in a R session


#@@@@@@@@@@@@@@@@#
##### INPUTS #####
#@@@@@@@@@@@@@@@@#
4 files(MANDATORY but extensible):
SEURAT_OUTPUT
STRELKA_SNV
STRELKA_INDEL
MUTECT_OUTPUT


#@@@@@@@@@@@@@@@@#
##### RUNNING ####
#@@@@@@@@@@@@@@@@#
To get usage run as : 
	bash Snipea.main.sh --help
or 
	Snipea.main.sh -h

Example with ALL the MANDATORY parameters:
bash ${DIR_INSTALL}/Snipea.main.sh \
--dirscript ${DIR_INSTALL} \
--seu  "my.seurat.vcf" \
--slksnv "my.strelka.passed.somatic.snvs.vcf" \
--slkindel "my.strelka.passed.somatic.indels.vcf" \
--mtcsnv "my.MuTect.filt.vcf"

optional parameters:
-h / --help : will show the usage message
--outprefix will assign a user's given name instead of automatically capturing the prefix from the Suerat's-Filtered-SNVs file.
--keepVennInputFiles : possible value: yes or no ; if yes, the input file for the Venn R script will not be deleted

#@@@@@@@@@@@@@@@@#
## INSTALLATION ##
#@@@@@@@@@@@@@@@@#
unzip the archive Snipea.Final.zip, and edit and modify the first line of the script named << Snipea.main.sh  >> with the appropriate path.
The path should be the full path of the installation directory. For example: export dir_script="${HOME}/tools/Snipea"


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
##### PREREQUISITES / INFOS  #####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
The assumption for this applications is that the 5 inputs files contains FILTERED variants, which will minimize the number of variants to process. If you'd like to change your own filtering parameters, please modify them in the "Snipea.SNV_filtering.sh"
See ## INPUTS ## section above for examples of expected files.

The PREFIX variable for the ouput filename of the merged file is based on the Seurat's input filename. 
The prefix is any string of character up to the first dot encounter into the filename.
for instance: if the filename is << my_prefix_from_filename.pass.seu.snvs.vcf >>, the prefix will be: << my_prefix_from_filename >>.

The bash scripts contains comments regarding the implementation and the choices made to merge the variant events from the 3 tools.

As Strelka calls the indels better, the prevalence for the INDELs has been assigned to Strelka. 
For the SNVs, the prevalence is Seurat, Strelak, Mutect, in that order.
Seurat gives two lines (or more) if a position has two (or more) different ALT. Strelka (and mutect) gives only one. 
If Strelka calls two alterations(ALT), as well as Seurat and Mutect; we create 1 line with defined Prevalence for SNVs.
If Strelka calls two ALT, BUT NOT Seurat and Mutect; we still create 1 line with defined Prevalence for SNVs, but the two ALT will still be kept.

If seurat at least, and all the other tools call that position, for each ALT, one line is written. Therefore, we will have two or more lines for a position with more than one ALT.
If seurat do not call that position, and still two ALT are found, the merger follows the same format as the tool has. For Strelka and mutect, they produce only one line. The ALT column will therefore have a list of ALT bases separated by a comma.

If Seurat is the only one having called the same position with two different alleles, we currently skip the positions;
for Strelka, if the ALT column contains a comma, meaning Two ALT alleles were called, we keep the ALT as is in the merge file, but we capture only the first allele to calculate AR and AD information.


Genotype Definition and thresholds:
Homozygote reference if AR<=0.200, heterozygote if 0.200<AR<0.900, and Homozygote ALT if AR>=0.900

We customized some FLAG and TAG; All of them have been added to the vcf Header.
The outputted vcf is VCF 4.1 compliant.


If any of the input vcf file has more than 10000 events within, an email will be sent to the user who started the script letting the user know that the merging might take a while. This threshold is hardcoded and can be modified by the owner of the script. (hard coded threshold in function <<countEventsPerFile>> )


GT has been captured based on the AR value. AR Thresholds have been hard coded to 0.200 and 0.900 in the functions "calculate_GT_using_AR_value" for both INDELS and SNVS;
look for lines: "local thld_min=0.200; local thld_max=0.900;" ; 

For Strelka, tiers 1 and 2 are available; To calculate AR, AD and get GT, only tier1 values have been used. tier1 values are slightly more stringent than tier2.

