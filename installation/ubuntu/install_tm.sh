#!/usr/bin/env bash

# Make sure to get exit status of last commant before the pipe
set -o pipefail


## install dependencies
echo "Installing system dependencies..."
source ${TM_HOME}/installation/ubuntu/install_dependencies.sh | tee ${TM_HOME}/.tm_install_deps.log

## Install all required R packages (or run Setup.R in RStudio)
echo "Installing required R packages"
${TM_HOME}/Setup.R | tee ${TM_HOME}/.tm_install_rsetup.log

## compile the parser needed to convert TissueAnalyzer outputs into csv
echo "Installing parser to read TissueAnlyzer outputs"
cd ${TM_HOME}/parser && make | tee ${TM_HOME}/.tm_install_parser.log
