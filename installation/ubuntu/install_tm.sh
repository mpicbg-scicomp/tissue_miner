#!/usr/bin/env bash

# Make sure to get exit status of last command before the pipe
set -o pipefail

## install dependencies
echo "Installing system dependencies..."
source ${TM_HOME}/installation/ubuntu/install_dependencies.sh | tee ${TM_HOME}/.tm_install_deps.log

if [ $? -eq 1 ]; then
  echo "Install dependencies FAILED, please check that > R3.2 is available for your OS"
  exit 1
fi

## Install all required R packages (or run Setup.R in RStudio)
echo "Installing required R packages"
sudo ${TM_HOME}/Setup.R | tee ${TM_HOME}/.tm_install_rsetup.log

if [ $? -eq 1 ]; then
  echo "The R setup FAILED, please check ${TM_HOME}/.tm_install_rsetup.log"
  exit 1
fi

## compile the parser needed to convert TissueAnalyzer outputs into csv
echo "Installing parser to read TissueAnlyzer outputs"
cd ${TM_HOME}/parser && make | tee ${TM_HOME}/.tm_install_parser.log

if [ $? -eq 1 ]; then
  echo "The parser compilation FAILED, please check ${TM_HOME}/.tm_install_parser.log"
  exit 1
fi


