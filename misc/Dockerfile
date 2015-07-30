
## see https://www.digitalocean.com/community/tutorials/docker-explained-using-dockerfiles-to-automate-building-of-images
## https://docs.docker.com/articles/dockerfile_best-practices/

## see https://registry.hub.docker.com/u/rocker/shiny/dockerfile/
FROM ubuntu

#RUN apt-get install -y git
RUN apt-get install -y curl

## install dependencies without cloning the repository. This allows for a more efficient use of the docker build cache
RUN /bin/bash -c "source <(curl https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/misc/install_dependencies.sh)"

## force image rebuild after changing sources by changing comment in line
RUN mkdir /tissue_miner && git clone https://github.com/mpicbg-scicomp/tissue_miner /tissue_miner ## test7

## http://stackoverflow.com/questions/20635472/using-the-run-instruction-in-a-dockerfile-with-source-does-not-work
#RUN /bin/bash -c "source /tissue_miner/misc/install_dependencies.sh"

## Install all required R packages
RUN  /bin/bash -c "export TM_HOME='/tissue_miner'; /tissue_miner/Setup.R | tee tm_setup.log"

## compile the parser needed to convert TissueAnalyzer outputs into csv
RUN cd /tissue_miner/parser &&  make

## also add .bash_profile to user home directory
COPY docker_bash_profile.sh /.bash_profile

## install virtual xserver because imageParser needs it
RUN apt-get install --assume-yes xvfb