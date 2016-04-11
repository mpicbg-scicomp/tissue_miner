#!/usr/bin/env bash

# Make sure to get exit status of last command before the pipe
set -o pipefail

## Ensure that tm installation directory is exported as TM_HOME
if [ -z "$TM_HOME" ]; then
    echo "TM_HOME is not defined" >&2
    exit 1
fi


## install dependencies
echo "Installing system dependencies..."
source ${TM_HOME}/installation/ubuntu/install_dependencies.sh | tee ${TM_HOME}/.tm_install_deps.log

if [ $? -eq 1 ]; then
  echo "Installation of dependencies failed" >&2
  exit 1
fi

## Install all required R packages (or run Setup.R in RStudio)
echo "Installing required R packages"
${TM_HOME}/Setup.R | tee ${TM_HOME}/.tm_install_rsetup.log

if [ $? -eq 1 ]; then
  echo "The R setup FAILED, please check ${TM_HOME}/.tm_install_rsetup.log" >&2
  exit 1
fi

## compile the parser needed to convert TissueAnalyzer outputs into csv
echo "Installing parser to read TissueAnlyzer outputs"
cd ${TM_HOME}/parser && make | tee ${TM_HOME}/.tm_install_parser.log

if [ $? -eq 1 ]; then
  echo "The parser compilation FAILED, please check ${TM_HOME}/.tm_install_parser.log" >&2
  exit 1
fi

## set up .tm_profile
echo "Set up .tm_profile"
cd $HOME

if [ ! -f .tm_profile ]; then 
  echo "## TissueMiner configuration:" >> .tm_profile
  echo "export TM_HOME=${TM_HOME}" >> .tm_profile
else
  grep -qx "export TM_HOME=${TM_HOME}" .tm_profile
  if [ $? -eq 1 ]; then 
    echo "## TissueMiner configuration:" >> .tm_profile
    echo "export TM_HOME=${TM_HOME}" >> .tm_profile
  fi
fi

diff .tm_profile ${TM_HOME}/installation/ubuntu/tm_profile.txt | grep "^>" | sed 's/^> //g' >> .tm_profile

## set up .bashrc
echo "Set up .bashrc"
if [ ! -f .bashrc ]; then
  echo "
# Source the .tm_profile
if [ -f ~/.tm_profile ]; then
   source ~/.tm_profile
fi
" > .bashrc
else
  grep -qw "if \[ -f ~/.tm_profile \]" .bashrc
  if [ $? -eq 1 ]; then 
    echo "
# Source the .tm_profile
if [ -f ~/.tm_profile ]; then
   source ~/.tm_profile
fi
" >> .bashrc
  fi
fi


echo "###\n###\nTM installation Done\n###\n###"