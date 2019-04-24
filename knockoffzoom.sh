#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Run KnockoffZoom on a toy dataset
#
# Authors: Matteo Sesia
# Date:    04/24/2019

########
# TODO #
########
# - Check R libraries
# - Do not load adjclust

######################
# Check dependencies #
######################
check_dependency () {
  CMD=$1
  if [ ! -x "$(command -v $CMD)" ]; then
    echo "Error: command $CMD not available"
    exit
  fi
}
DEPENDENCY_LIST=("plink" "plink2" "fastphase" "datamash" "R")
for DEPENDENCY in "${DEPENDENCY_LIST[@]}"; do  
  check_dependency $DEPENDENCY
done

#########
# Setup #
#########

# Setup spinner for long jobs
source "knockoffzoom/utils/spinner.sh"

####################
# Run KnockoffZoom #
####################

# Print banner
cat "knockoffzoom/banner.txt"

# Log file
LOG_FILE="knockoffzoom.log"
rm -r $LOG_FILE
echo "Output logged in: "$LOG_FILE

# Enter the directory where the scripts are stored
cd knockoffzoom

# Module 1: fit the HMM
start_spinner 'Running module 1...'
sleep 1
#./module_1_hmm.sh &> "../"$LOG_FILE
stop_spinner $?

# Module 2: partition the genome into LD blocks
start_spinner 'Running module 2...'
sleep 1
#./module_2_partition.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 3: generate the knockoffs
start_spinner 'Running module 3...'
./module_3_knockoffs.sh &>> "../"$LOG_FILE
stop_spinner $?

exit

# Module 4: compute the test statistics
start_spinner 'Running module 4...'
./module_4_statistics.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 5: report significant findings
start_spinner 'Running module 4...'
./module_5_discover.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 6: plot the results (TODO)
