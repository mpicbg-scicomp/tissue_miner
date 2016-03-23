########################################################################################################################
## Install requirements for TissueMiner

## (optionally) install R 3.2.1 (see http://www.thertrader.com/2014/09/22/installing-rrstudio-on-ubuntu-14-04/)
if [ -z "$(which R)" ]; then
echo "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -sc)/" | sudo tee --append /etc/apt/sources.list
#sudo  cat /etc/apt/sources.list

# work around permission issues
# see http://stackoverflow.com/questions/10255082/installing-r-from-cran-ubuntu-repository-no-public-key-error
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
sudo apt-get install -y r-base
fi

## enforce at least R 3.2
if [ -z "$(R --version | grep '3.2')" ]; then
    echo "R is too old. At least v3.2 is required to run TissueMiner" >&2
	echo "Please, try to update it by running: sudo apt-get update && sudo apt-get install -y r-base" >&2
    return 1
fi

## install snakemake
sudo apt-get install -y python3-setuptools
easy_install3 snakemake

## install curl
sudo apt-get --no-install-recommends install -y libcurl4-gnutls-dev

## install lixbml
sudo apt-get --no-install-recommends install -y libxml2-dev

## install imagemagick needed for image rotation
sudo apt-get --no-install-recommends install -y bc imagemagick

## install openssl for git
sudo apt-get --no-install-recommends install -y libssl-dev

## install rgeos needed for polygon operations
sudo apt-get --no-install-recommends install -y libgeos-dev

## install sqlite
sudo apt-get --no-install-recommends install -y sqlite3 libsqlite3-dev

## install avconv needed movie rendering
sudo apt-get --no-install-recommends install -y libav-tools

## install graphviz needed for snakemake graph rendering
sudo apt-get --no-install-recommends install -y graphviz

## needed for svg export via sq
sudo apt-get --no-install-recommends install -y libcairo2-dev


## sem is disabled for now, since it's not essential to have it.
## sem for image conversion (see http://askubuntu.com/questions/12764/where-do-i-get-a-package-for-gnu-parallel)
#sudo apt-get install -y parallel
#sudo apt-get install -y wget
#wget http://ftp.gnu.org/gnu/parallel/parallel-20140422.tar.bz2
#tar -xvjf parallel*
#cd parallel*
##less README
#./configure
#make clean
#make


#sudo make install
#./src/sem --bibtex
# note added to PATH via .bash_profile in docker image

#### install image parser dependencies

# netcds http://askubuntu.com/questions/79418/installing-netcdf
#sudo apt-get install -y libnetcdf-dev
# see http://askubuntu.com/questions/382444/qt-headers-and-libraries-not-found
sudo apt-get install -y qt4-dev-tools libqt4-core
#sudo apt-get install libqt4-dev libqt4-opengl-dev libqtwebkit-dev qt4-linguist-tools qt4-qmake

#
### git should be already present for cloning TM from git hub, except if TM was just manually downloaded
sudo apt-get --no-install-recommends install -y git
