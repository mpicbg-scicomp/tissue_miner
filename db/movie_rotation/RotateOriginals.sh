#!/bin/bash

if [ $# -ne 1 ]; then echo "Usage: `basename $0` <movie db directory>"; exit; fi

# http://stackoverflow.com/questions/3915040/bash-fish-command-to-print-absolute-path-to-a-file
export movieDbDir=$(readlink -e $1)
    #export movieDbDir=$(pwd)

# Half nb cpu on Linux machines
NB_CPU=$(echo $(grep -c processor /proc/cpuinfo)/2 | bc)
#NB_CPU="4"

# half nb cpu on mac
#NB_CPU=$(echo $(/usr/sbin/system_profiler -detailLevel full SPHardwareDataType | awk '/Total Number of Cores/{print $5}')/2 | bc)

### fixme not generic
### cluster support
export PATH=/projects/project-raphael/parallel-20140422/src:$PATH
### bioinfo support
export PATH=/home/brandl/bin/parallel-20140422/src:$PATH

#export trafoFile="$movieDbDir/transformation.txt"
export trafoFile="$movieDbDir/Segmentation/transformation.txt"

if [ ! -f "$trafoFile" ]; then
    echo "no transformation.txt, skipping png rotation"; exit
fi

for originalPng in $(find $movieDbDir/Segmentation -name "original.png" | sort); do
    # DEBUG originalPng=/lustre/projects/project-raphael/movie_dbs/testing/WT_25deg_111103/Segmentation/frame0204/original.png

    trafoPngOutput=${originalPng%%.png}_trafo.png

    echo "rotating $originalPng to $trafoPngOutput ..."
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