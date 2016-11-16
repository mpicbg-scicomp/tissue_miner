TM installation directly on your operating system (not recommended)
================


### Ubuntu

TissueMiner ships with a command-line Ubuntu-Installer

* To install TissueMiner on Ubuntu simply copy-paste the commands in the box below into a [terminal](https://help.ubuntu.com/community/UsingTheTerminal) to perform the following steps:
    + Define the directory where to install TissueMiner
    + Download TissueMiner repository (~120MB)
    + Install TM (~30 min due to the compilation of the source code of R packages)
```
# Step1
export TM_HOME="${HOME}/tissue_miner"
sudo apt-get update && sudo apt-get install git

```
```
# Step 2
git clone https://github.com/mpicbg-scicomp/tissue_miner.git ${TM_HOME}
${TM_HOME}/installation/ubuntu/install_tm.sh
source ${HOME}/.bashrc
```

* To update TissueMiner on Ubuntu, click [here](../faq.md#how-to-update-my-tissueminer-installation)


### MacOS

