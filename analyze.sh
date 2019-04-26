#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Run KnockoffZoom on a toy dataset
#
# Authors: Matteo Sesia
# Date:    04/24/2019

#########
# Setup #
#########

# Print header
printf "KnockoffZoom v0.1 (24 Apr 2019) \n"
printf "https://bitbucket.org/msesia/knockoffzoom \n"
printf "(C) 2019 Matteo Sesia (Stanford University)   GNU General Public License v3 \n\n"

# Setup spinner for long jobs
source "misc/spinner.sh"

# Log file
LOG_FILE="knockoffzoom.log"
rm -f $LOG_FILE
touch $LOG_FILE
echo "Log file: "$LOG_FILE

# Enter the directory where the scripts are stored
cd knockoffzoom

######################
# Check dependencies #
######################

printf "\nSetup\n"
# System dependencies
ERROR=0
check_dependency () {
  CMD=$1
  if [ ! -x "$(command -v $CMD)" ]; then
    echo -e "Error: command $CMD not available"
    ERROR=1
  fi
}
DEPENDENCY_LIST=("plink" "plink2" "fastphase" "datamash" "awk" "R")
start_spinner " - Checking system dependencies..."
for DEPENDENCY in "${DEPENDENCY_LIST[@]}"; do  
  check_dependency $DEPENDENCY &>> "../"$LOG_FILE
done
stop_spinner $ERROR

# R libraries
start_spinner " - Checking R library dependencies..."
Rscript --vanilla "utils/check_packages.R" &>> "../"$LOG_FILE
stop_spinner $?

####################
# Run KnockoffZoom #
####################

printf "\nData analysis\n"

# Module 1: fit the HMM
start_spinner ' - Running module 1...'
./module_1_hmm.sh &> "../"$LOG_FILE
stop_spinner $?

# Module 2: partition the genome into LD blocks
start_spinner ' - Running module 2...'
./module_2_partition.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 3: generate the knockoffs
start_spinner ' - Running module 3...'
./module_3_knockoffs.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 4: compute the test statistics
start_spinner ' - Running module 4...'
./module_4_statistics.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 5: report significant findings
start_spinner ' - Running module 5...'
./module_5_discover.sh &>> "../"$LOG_FILE
stop_spinner $?

#####################
# Summarize results #
#####################

printf "\nResults written in 'results/'\n"
