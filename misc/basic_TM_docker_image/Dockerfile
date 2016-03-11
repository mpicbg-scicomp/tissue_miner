
## see https://www.digitalocean.com/community/tutorials/docker-explained-using-dockerfiles-to-automate-building-of-images
## https://docs.docker.com/articles/dockerfile_best-practices/

## see https://registry.hub.docker.com/u/rocker/shiny/dockerfile/
FROM ubuntu

MAINTAINER "Raphael Etournay and Holger Brandl" brandl@mpi-cbg.de

#RUN apt-get install -y git
RUN apt-get install -y curl

## install openssl for git
#RUN apt-get install -y libssl-dev
RUN apt-get install -y git


## install dependencies without cloning the repository. This allows for a more efficient use of the docker build cache
RUN /bin/bash -c "source <(curl https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/misc/basic_TM_docker_image/install_dependencies.sh)"  ## bump1

## force image rebuild after changing sources by changing comment in line
RUN mkdir /tissue_miner && git clone https://github.com/mpicbg-scicomp/tissue_miner /tissue_miner ## bump1

### http://stackoverflow.com/questions/20635472/using-the-run-instruction-in-a-dockerfile-with-source-does-not-work
##RUN /bin/bash -c "export TM_HOME='/tissue_miner'; source /tissue_miner/install_tm.sh"

## Install all required R packages
RUN  /bin/bash -c "export TM_HOME='/tissue_miner'; /tissue_miner/Setup.R | tee ${TM_HOME}/.tm_install_rsetup.log"

## compile the parser needed to convert TissueAnalyzer outputs into csv
RUN cd /tissue_miner/parser &&  make


## also add .bash_profile to user home directory
COPY docker_bash_profile.sh /.bash_profile

## install virtual xserver because imageParser needs it
# run without x https://linuxmeerkat.wordpress.com/2014/10/17/running-a-gui-application-in-a-docker-container/
#RUN apt-get install -y xvfb


COPY docker_entrypoint.sh /
ENTRYPOINT ["/docker_entrypoint.sh"]
#CMD ["sm", "-n"]
