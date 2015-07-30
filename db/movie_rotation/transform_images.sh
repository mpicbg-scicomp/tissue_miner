#!/bin/sh

if [ $# -ne 3 ]; then echo "Usage: `basename $0` <inFilePath> <transfoPath> <outFilePath>" >&2 ; exit; fi


## ROTATE (AND FLIP) IMAGES BY READING PARAMETERS FROM transformation.txt
## REQUIRES imagemagick to be installed

export PATH=$PATH:/projects/project-raphael

trafoImage(){
	if [ $3 -eq 0 ]
	then
		echo "Rotate $1 --> $4"		
		convert $1 -background "black" -rotate "$2" -type truecolor -colorspace RGB -define png:compression-level=9 $4
	else
		echo "Rotate and flip $1 --> $4"
		convert $1 -background "black" -rotate "$2" -flip -type truecolor -colorspace RGB -define png:compression-level=9 $4
	fi
}


## Input variables
inFilePath=$1
transfoPath=$2
outFilePath=$3

## DEBUG
#inFilePath=/home/etournay/RawData/original.png
#transfoPath=/home/etournay/RawData/transformation.txt
#outFilePath=/home/etournay/RawData/toto.png


if [ ! -f $inFilePath ] || [ ! -f $transfoPath ] || [ ! -d $(dirname $outFilePath) ]
then
	echo "Missing input, aborting..."; #exit 0
else
	anglerad=$(cat $transfoPath | awk 'NR>1{print $1}')
	pi=$(echo "scale=10; 4*a(1)" | bc -l)
	angledeg=$(echo "scale=3; $anglerad*180/$pi" | bc -l)
	IsVerticalFlip=$(cat $transfoPath | awk 'NR>1{print $6}')
	#echo "anglerad=$anglerad; pi=$pi; angledeg=$angledeg; IsVerticalFlip=$IsVerticalFlip"
	trafoImage $inFilePath $angledeg $IsVerticalFlip $outFilePath
fi

exit 0


