#! /bin/bash

## Avconv installation 
hash avconv &> /dev/null
if [ $? -eq 1 ]; then
    brew install faac lame speex libogg libvorbis opencore-amr freetype theora x264 yasm
    brew install libav
fi

## R installation 
hash R &> /dev/null
if [ $? -eq 1 ]; then
    brew tap homebrew/science
    brew install r 
fi

## libgeos_c library installation
# http://collectivegeo.readthedocs.org/en/latest/installation.html
brew install geos

## Configure R by running setup.R
${TM_HOME}/Setup.R


# Rstudio installation
brew install Caskroom/cask/rstudio # OK, installed in $HOME/Applications/

