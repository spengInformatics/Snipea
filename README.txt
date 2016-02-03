#!/bin/bash
#####################################################################
##### Author: Sen Peng #####
##### Date: December 2014#####
#######################################
#@@@@@@@@@@@@@@@@#
## INTRODUCTION ##
#@@@@@@@@@@@@@@@@#

Snipea (SNv Integration, Prioritization, Ensemble, and Annotation) is a tool to further identify, prioritize and annotate somatic SNV from patient tumor (Other Biological samples). Our tool will apply an ensemble-based approach to improve the accuracy and annotation of variant callers. “Snipea” calculates the authenticity and statistical significance that an SNV is causative for a query disease and hence provides a means of prioritizing candidate SNVs. The ranking algorithm will take into account the callers_agreement, clinical_impact (protein-altering aberrations such as non-synonymous SNV, Stop code gained/loss or Frame shift mutation), quality score from each tool and public_database_status.

#@@@@@@@@@@@@@@@@#
## REQUIREMENTS ##
#@@@@@@@@@@@@@@@@#
This software requires the following tools/software to be installed


linux commands: cat, case, cd , cut, sed, grep, awk, sort, uniq, getopt, wc, mv, cp, head , pwd, echo, for while, if, bash, read, mkdir, wait, do, rm, local.
bedtools version 2.17 or up
R packages for version 2.15.2 of R: 
		VennDiagram 1.6.5 or up
		Vennerable 2.2 or up
		gplots 2.11.0 or up

#@@@@@@@@@@@@@@@@#
##### RUNNING ####
#@@@@@@@@@@@@@@@@#
To get usage run as : 
	bash Snipea.main.sh --help
or 
	Snipea.main.sh -h





