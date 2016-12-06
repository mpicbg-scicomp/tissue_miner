#! /bin/bash
echo "Fetching TM docker-image version on github..."
TM_GITVERSION=$(curl -s https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/dockerimg_version.txt)

# Test if TM version isn't empty
if [ -z $TM_GITVERSION ]; then
  echo "Fetching TM docker-image version on github failed, skipping..." 
  exit 1
fi

# Parse TM version from git-hub
TM_GITMAJOR=$(echo $TM_GITVERSION | awk -F "." '{print $1}' | bc)
TM_GITMINOR=$(echo $TM_GITVERSION | awk -F "." '{print $2}' | bc)
TM_GITREV=$(echo $TM_GITVERSION | awk -F "." '{print $3}' | bc)

# Test if one variable is empty 
if [ -z $TM_GITMAJOR ] || [ -z $TM_GITMINOR ] || [ -z $TM_GITREV ]; then 
echo "Abnormal TM docker-image version on github, skipping..." 
  exit 1
fi

# Parse local TM version 
TM_LOCVERSION=$(cat $TM_HOME/dockerimg_version.txt)
TM_LOCMAJOR=$(echo $TM_LOCVERSION | awk -F "." '{print $1}')
TM_LOCMINOR=$(echo $TM_LOCVERSION | awk -F "." '{print $2}' | bc)
TM_LOCREV=$(echo $TM_LOCVERSION | awk -F "." '{print $3}' | bc)

# Helper function to display the update message
function display_message {
RED="\033[1;31m"
NO_COLOUR="\033[0m"
echo -e "${RED}Your current TM_${TM_LOCVERSION} is OUTDATED !\n\
Please UPDATE to TM_${TM_GITVERSION} using: \
'docker pull etournay/tissue_miner${NO_COLOR}'"
}

# Check if TM needs to be updated
if [ $TM_GITMAJOR -gt $TM_LOCMAJOR ]; then
	display_message
	else if [ $TM_GITMAJOR -eq $TM_LOCMAJOR ] && [ $TM_GITMINOR -gt $TM_LOCMINOR ]; then
		display_message
		else if [ $TM_GITMINOR -eq $TM_LOCMINOR ] && [ $TM_GITREV -gt $TM_LOCREV ]; then
			display_message
		fi
	fi
fi




