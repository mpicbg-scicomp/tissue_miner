#!/usr/bin/env bash


## install dependencies
source ${TM_HOME}/misc/install_dependencies.sh

## Install all required R packages (or run Setup.R in RStudio)
${TM_HOME}/Setup.R | tee tm_setup.log

## compile the parser needed to convert TissueAnalyzer outputs into csv
cd ${TM_HOME}/parser && make

## adjust your path to include all tools
export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
export PATH=${TM_HOME}/parser:$PATH