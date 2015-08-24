#!/bin/bash

if [ $# -ne 1 ]; then echo "Usage: `basename $0` <Movie Data Directory>" >&2 ; exit; fi

export movieDbDir=$1
#export movieDbDir=$(pwd)

if [ ! -d "$movieDbDir" ]; then  echo "directory $movieDbDir does not exist or is not a directory" >&2  exit; fi


movieName=$(basename $movieDbDir)
segDataDir=$movieDbDir/Segmentation

## create the config file for the image parser
#imFile=$(find $(dirname $segDataDir) -name "$movieName*.[p,t][n,i][f,g]" | head -n1)
imFile=$(find $(dirname $segDataDir) -name "tracked_cells_resized.tif" | head -n1)
imWidth=$(identify $imFile | awk '{print $3}' | awk 'BEGIN{FS="x"}{print$1}')
imHeight=$(identify $imFile| awk '{print $3}' | awk 'BEGIN{FS="x"}{print$2}')
#numFrames=$(find $segDataDir | awk '/tracked_cells_resized.tif/{print}'| wc -l)
numFrames=$(find $segDataDir -name "tracked_cells_resized.tif" | wc -l)

echo "Copy $imWidth $imHeight $numFrames 16 0.0834 1 0.208 into movieInformation.dat"
echo "$imWidth $imHeight $numFrames 16 0.0834 1 0.208" > $segDataDir/movieInformation.dat
cat $segDataDir/movieInformation.dat

exit 0
