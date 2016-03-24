#! /bin/bash

# Set colors
BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m'

# Display a message for the user
echo ""
echo -e ${BLUE}"Welcome to the TissueMiner installation procedure"${NC}
echo ""
echo "This will install all required dependencies in your account without affecting already existing system tools such as Virtualbox:"
echo "  - Homebrew"
echo "  - VirtualBox"
echo "  - Docker Machine"
echo "  - Docker"
echo "  - Dockerized TissueMiner"
echo ""
echo "... and it may take a while depending on your OS configuration and Internet connection"
echo ""
echo "You may be asked to type in your user password during the installation process"
echo ""

# Get computer hw info and set paramters for the docker virtual machine
SYST_NB_CPU=$(sysctl -n hw.ncpu)
SYST_MEM_SIZE=$(sysctl -n hw.memsize)

VM_NAME="tm-vm"
DOCKER_IMAGE="etournay/tissue_miner"
# DOCKER_IMAGE="ubuntu" # for DEBUG

VM_NB_CPU=$(echo $(( SYST_NB_CPU/2 )) | bc)
VM_MEM_SIZE=$(echo $(( SYST_MEM_SIZE/(2*1024*1024) )) | bc)


# Install Homebrew package manager and update DB
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
echo ""
echo "Update brew database, it may take a while..."
echo ""
brew update


# Docker installation including virtualbox, docker-machine and docker
brew cask install dockertoolbox # CAUTION requires to give a password

# Check Docker installation
if [ -z "$(which docker-machine)" ]; then
  echo "Docker Machine is not found. Please check if any error occurred earlier during installation."
  read -r -p "Press Enter to quit the installer..." key
  exit 1
fi

if [ -z "$(which VBoxManage)" ]; then
  echo "VirtualBox is not found. Please check if any error occurred earlier during installation."
  read -r -p "Press Enter to quit the installer..." key
  exit 1
fi

# Check if a VM called $VM_NAME exists
VBoxManage list vms | grep \""${VM_NAME}"\" &> /dev/null
VM_EXISTS_CODE=$?

# If the VM doesn't exit then create it
if [ $VM_EXISTS_CODE -eq 1 ]; then
  docker-machine rm -f "${VM_NAME}" &> /dev/null
  rm -rf ~/.docker/machine/machines/"${VM_NAME}"
  # That creates a VM using 50% of computer CPUs and 50% of computer memory.
  docker-machine create --virtualbox-cpu-count "${VM_NB_CPU}" --virtualbox-memory "${VM_MEM_SIZE}" --virtualbox-disk-size 204800 -d virtualbox "${VM_NAME}"
fi

VM_STATUS="$(docker-machine status ${VM_NAME} 2>&1)"
if [ "${VM_STATUS}" != "Running" ]; then
  docker-machine start "${VM_NAME}"
  yes | docker-machine regenerate-certs "${VM_NAME}"
fi

# Load the VM environment into the current shell for enabling docker
eval "$(docker-machine env ${VM_NAME})"

# Display the docker 
clear
cat << EOF


                        ##         .
                  ## ## ##        ==
               ## ## ## ## ##    ===
           /"""""""""""""""""\___/ ===
      ~~~ {~~ ~~~~ ~~~ ~~~~ ~~~ ~ /  ===- ~~~
           \______ o           __/
             \    \         __/
              \____\_______/


EOF
echo -e "${BLUE}docker${NC} is configured to use the ${GREEN}${VM_NAME}${NC} machine with IP ${GREEN}$(docker-machine ip ${VM_NAME})${NC}"

# Get the VM state (active or not)
#VM_STATE=$(docker-machine ls | awk '/'${VM_NAME}'/{print $2}')

# Pull dockerized TissueMiner according to VM state
  #https://github.com/docker/docker/issues/1888

#if [ "${VM_STATE}" != "*" ]; then
#  docker-machine ssh "${VM_NAME}" docker pull ${DOCKER_IMAGE}
#else
#  docker pull ${DOCKER_IMAGE}
#fi

# Check if the TissueMiner image is present
#if [ "${VM_STATE}" != "*" ]; then
#  echo ""
#  echo "Your docker environment is not accessible without ssh"
#  echo ""
#  echo -e "Please, make sure that the tissue_miner image is present with the command:"
#  echo "docker-machine ssh ${VM_NAME} docker images"
#  echo ""
#else
#  IMG=$(docker images | awk '/'${DOCKER_IMAGE}'/{print $1}')
#  if [ "${IMG}" != "${DOCKER_IMAGE}" ];then
#    echo ""
#    echo "${DOCKER_IMAGE} not found. The TissueMiner installation failed !"
#    echo "Please, try again..."
#    echo ""
#  else
#    echo ""
#    echo "A dockerized TissueMiner ${DOCKER_IMAGE} has been successfully installed and configured on your computer."
#    echo ""
#  fi
#fi

read -r -p "Press Enter to quit the installer..." key


