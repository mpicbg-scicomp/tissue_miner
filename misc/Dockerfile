
## see https://www.digitalocean.com/community/tutorials/docker-explained-using-dockerfiles-to-automate-building-of-images
## https://docs.docker.com/articles/dockerfile_best-practices/

## see https://registry.hub.docker.com/u/rocker/shiny/dockerfile/
FROM ubuntu

RUN apt-get install -y git

## force image rebuild after changing sources by changing comment in line
RUN mkdir /tissue_miner && git clone https://github.com/mpicbg-scicomp/tissue_miner /tissue_miner/ ## test3

## http://stackoverflow.com/questions/20635472/using-the-run-instruction-in-a-dockerfile-with-source-does-not-work
RUN /bin/bash -c "source /tissue_miner/misc/install_dependencies.sh"

#RUN ${TM_HOME}/tissue_miner/Setup.R | tee tm_setup.log

## compile the parser needed to convert TissueAnalyzer outputs into csv
#RUN cd ${TM_HOME}/parser && make clean && make

## also add .bash_profile
COPY docker_bash_profile.sh /.bash_profile