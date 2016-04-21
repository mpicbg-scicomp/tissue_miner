#!/bin/bash


#while getopts t flag; do
#  case $flag in
#    t)
#      echo "Preparing movie db";
#      useRotation=TRUE
#      ;;
#    ?)
#      exit;
#      ;;
#  esac
#done
#shift $(( OPTIND - 1 ));


if [ $# -ne 1 ]; then echo "Usage: `basename $0` <Movie Data Directory>" >&2 ; exit; fi

export movieDbDir=$1

#DEBUG movieDbDir=/mnt/project_raphael@filesServer/movieSegmentation/WT_25deg_111102
#DEBUG movieDbDir=/mnt/project_raphael@fileserver/movieSegmentation/120531_Debug_Tissue_Sample
#DEBUG movieDbDir=/mnt/project_raphael@fileserver/movieSegmentation/PA_Sample_NoCorrection


if [ ! -d "$movieDbDir" ]; then  echo "directory $movieDbDir does not exist" >&2  exit; fi


#echo "done using rotation=$useRotation"
#exit

movieName=$(basename $movieDbDir)
#cd $movieName


segDataDir=$movieDbDir/Segmentation


########################################################################################################################
## create original.png from tif or png time-lapse images as this is required for the parser, and not generated in TA anymore
for originalIm in $(find $segDataDir -name "$movieName*.[p,t][n,i][f,g]" | sort); do
        # DEBUG originalPng=/lustre/projects/project-raphael/movie_dbs/testing/WT_25deg_111103/image_data/mutant/tag/segmentationData/frame0204/original.png
		file=$(basename $originalIm)		
		ext="${file##*.}"
        pngOutput=$(dirname $originalIm)/${file%.$ext}/"original".png
        echo "converting $originalIm -> $pngOutput ..."
	#sem -j4 convert $originalIm -background "black" -type truecolor -colorspace RGB -define png:compression-level=9 $pngOutput
	convert $originalIm -background "black" -type truecolor -colorspace RGB -define png:compression-level=9 $pngOutput
    done
#sem --wait


### create fake config file for parser NO ANYMORE NEEDED !!
#echo "Touch cumultimesec_for_PA.txt"Ω
#touch  $segDataDir/cumultimesec_for_PA.txt
#cp /projects/project_Raphael/A_Holger_Rapha_project/Rcode/commons/parser/movieInformation.dat $segDataDir/..

## create the config file for the image parser
#imFile=$(find $(dirname $segDataDir) -name "$movieName*.[p,t][n,i][f,g]" | head -n1)
imFile=$(find $(dirname $segDataDir) -name "tracked_cells_resized.tif" | head -n1)
imWidth=$(identify $imFile | awk '{print $3}' | awk 'BEGIN{FS="x"}{print$1}')
imHeight=$(identify $imFile| awk '{print $3}' | awk 'BEGIN{FS="x"}{print$2}')
numFrames=$(find $segDataDir -name "tracked_cells_resized.tif" |wc -l)

echo "Copy $imWidth $imHeight $numFrames 16 0.0834 1 0.208 into movieInformation.dat"
echo "$imWidth $imHeight $numFrames 16 0.0834 1 0.208" > $segDataDir/movieInformation.dat
cat $segDataDir/movieInformation.dat




exit 0
