#!/bin/bash
#############################################
# execute script in current directory
#$ -cwd
#$ -e HCCp1_Step3_4_DC_subset_individual_sorting.STDOUT.log
#$ -o HCCp1_Step3_4_DC_subset_individual_sorting.STDERR.log
# shell for qsub to use:
#$ -S /bin/bash
# name for the job; used by qstat
#$ -N HCCp1_Step3_4_DC_subset_individual_sorting

# only run on a node with given minimum memory free

## #$ -l mem_free=33G

# kill this job if its memory usage exceeds this limit

## #$ -l h_vmem=50G

######################################################
# specify all parameters to this qsub script in double-quotes
# Example
# qsub HCCp1_Step3_4_DC_subset_individual_sorting "param1 param2 param3 .."
######################################################

######################################################
# Append environmental varibles required to run program
# Examples:
# export PATH=~/koelschnj/HCC_paper1:${PATH}
# source /usr/local/bar/enable_bar.sh
######################################################
source /usr/global/R-4.1.2/enable_R.sh

# Variable used to measure time 
START=$(date +%s)
NSTART=$(date +%s.%N)

Rscript ~/HCC_paper1/Analysis_Scripts/Rscript_HCCp1_Step3_4_DC_subset_individual_sorting.R


# Some information about the node this job ran on:

echo ""
echo "memory usage: "
/usr/bin/free

echo Time is $(date)
echo Directory is $(pwd)
echo Path is $PATH
echo Ld Library Path is $LD_LIBRARY_PATH

echo ""
echo "R version: "

which R
R --version


