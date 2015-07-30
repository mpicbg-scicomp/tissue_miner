########################################################################################################################
## Install requirements for TissueMiner

## install R 3.2.1 (see http://www.thertrader.com/2014/09/22/installing-rrstudio-on-ubuntu-14-04/)
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
sudo  cat /etc/apt/sources.list

# work around permission issues
# see http://stackoverflow.com/questions/10255082/installing-r-from-cran-ubuntu-repository-no-public-key-error
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
sudo apt-get install --assume-yes r-base

## install snakemake
sudo apt-get install --assume-yes python3-setuptools
easy_install3 snakemake

## install curl
sudo apt-get install --assume-yes libcurl4-gnutls-dev

## install lixbml
sudo apt-get install --assume-yes libxml2-dev

## install imagemagick needed for image rotation
sudo apt-get install --assume-yes bc imagemagick

## install openssl for git
sudo apt-get install --assume-yes libssl-dev

## install rgeos needed for polygon operations
sudo apt-get install --assume-yes  libgeos-dev

## install sqlite
sudo apt-get install --assume-yes sqlite3 libsqlite3-dev

## install avconv needed movie rendering
sudo apt-get install --assume-yes libav-tools

## sem for image conversion (see http://askubuntu.com/questions/12764/where-do-i-get-a-package-for-gnu-parallel)
#sudo apt-get install --assume-yes parallel
sudo apt-get install --assume-yes wget
wget http://ftp.gnu.org/gnu/parallel/parallel-20140422.tar.bz2
tar -xvjf parallel*
cd parallel*
#less README
./configure
make clean
make
#sudo make install
#./src/sem --bibtex
# note added to PATH via .bash_profile in docker image

#### install image parser dependencies

# netcds http://askubuntu.com/questions/79418/installing-netcdf
sudo apt-get install --assume-yes libnetcdf-dev
# see http://askubuntu.com/questions/382444/qt-headers-and-libraries-not-found
sudo apt-get install --assume-yes qt4-dev-tools libqt4-core
#sudo apt-get install libqt4-dev libqt4-opengl-dev libqtwebkit-dev qt4-linguist-tools qt4-qmake


## git
sudo apt-get install --assume-yes git


########################################################################################################################
## Install TissueMiner

# todo adjust to local setup
export TM_HOME="/tissue_miner/"

mkdir ${TM_HOME}
## fixme use github url here
git clone ssh://brandl@fileserver.mpi-cbg.de/projects/project_Raphael/scripts/.flywing_base_git_origin ${TM_HOME}

## update existing local copy
#cd  ${TM_HOME}; git pull origin

## Run setup routine to install all require packages
${TM_HOME}/Setup.R 2>&1 | tee tm_setup.log


## todo make sure that log file does not contain any errors

## compile the image parser
cd $TM_HOME/parser
make clean
make
