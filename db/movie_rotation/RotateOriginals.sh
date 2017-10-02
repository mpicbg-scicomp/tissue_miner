#!/bin/bash

if [ $# -ne 1 ]; then echo "Usage: `basename $0` <movie db directory>"; exit; fi

# http://stackoverflow.com/questions/3915040/bash-fish-command-to-print-absolute-path-to-a-file
export movieDbDir=$(readlink -e $1)

# Define number of CPUs on Linux or MacOSX OS
unamestr=$(uname)
if [[ "$unamestr" == 'Linux' ]]; then
   NB_CPU=$(echo $(grep -c processor /proc/cpuinfo)/2 | bc)
elif [[ "$unamestr" == 'Darwin' ]]; then
   NB_CPU=$(echo $(/usr/sbin/system_profiler -detailLevel full SPHardwareDataType | awk '/Total Number of Cores/{print $5}')/2 | bc)
fi

# Define path to transformation.txt file and test if it exists
export trafoFile="$movieDbDir/Segmentation/transformation.txt"

if [ ! -f "$trafoFile" ]; then
    echo "no transformation.txt, skipping png rotation"; exit
fi


# Rotate and/or flip-flop images using ImageMagick
for originalPng in $(find $movieDbDir/Segmentation -name "original.png" | sort); do

    trafoPngOutput=${originalPng%%.png}_trafo.png

    #echo "rotating $originalPng to $trafoPngOutput ..."
    #sem --no-notice -j6  $TM_HOME/db/movie_rotation/transform_images.sh $originalPng $trafoFile $trafoPngOutput

    if [ -z "$(which sem)" ]; then
        $TM_HOME/db/movie_rotation/transform_images.sh $originalPng $trafoFile $trafoPngOutput
    else
        sem -j$NB_CPU $TM_HOME/db/movie_rotation/transform_images.sh $originalPng $trafoFile $trafoPngOutput
    fi
done

if [ ! -z "$(which sem)" ]; then
    sem --wait
fi
