#!/usr/bin/env bash


## install dependencies
source ${TM_HOME}/misc/install_dependencies.sh

## Install all required R packages (or run Setup.R in RStudio)
${TM_HOME}/Setup.R | tee tm_setup.log

## compile the parser needed to convert TissueAnalyzer outputs into csv
cd ${TM_HOME}/parser && make
