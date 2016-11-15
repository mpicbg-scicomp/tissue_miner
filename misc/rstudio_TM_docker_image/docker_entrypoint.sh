#!/bin/bash
set -e

eval "source /home/rstudio/.bash_profile; $@"

$TM_HOME/misc/./version.sh
