## Rstudio server setup on Ubuntu
# Inspired from https://github.com/rocker-org/rocker
# 486a544e072d


## see https://www.digitalocean.com/community/tutorials/docker-explained-using-dockerfiles-to-automate-building-of-images
# https://docs.docker.com/articles/dockerfile_best-practices/

# see https://registry.hub.docker.com/u/rocker/shiny/dockerfile/


FROM ubuntu:14.04

MAINTAINER "Raphael Etournay and Holger Brandl" brandl@mpi-cbg.de


## security updates for ubuntu and install curl and git
RUN apt-get update && apt-get upgrade -y && apt-get install -y curl nano wget htop

## install dependencies without cloning the repository. This allows for a more efficient use of the docker build cache
RUN /bin/bash -c "source <(curl https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/installation/ubuntu/install_dependencies.sh | sed 's/sudo//g')"  ## bump2


## force image rebuild after changing sources by changing comment in line
RUN mkdir /tissue_miner && git clone https://github.com/mpicbg-scicomp/tissue_miner /tissue_miner ## bump

### http://stackoverflow.com/questions/20635472/using-the-run-instruction-in-a-dockerfile-with-source-does-not-work
##RUN /bin/bash -c "export TM_HOME='/tissue_miner'; source /tissue_miner/installation/ubuntu/install_tm.sh"

## install virtual xserver because imageParser needs it
# run without x https://linuxmeerkat.wordpress.com/2014/10/17/running-a-gui-application-in-a-docker-container/
#RUN apt-get install -y xvfb

## compile the parser needed to convert TissueAnalyzer outputs into csv
RUN cd /tissue_miner/parser &&  make



## Set a default user. Available via runtime flag `--user rstudio` 
## Add user to 'rstudio' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (for rstudio or linked volumes to work properly). 

# Step 3
RUN useradd rstudio \
	&& mkdir /home/rstudio \
	&& chown rstudio:rstudio /home/rstudio \
	&& addgroup rstudio staff


## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
# Step 5
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

# Step 6
ENV LC_ALL en_US.UTF-8
# Step 7
ENV LANG en_US.UTF-8


## Add RStudio binaries to PATH
# Step 9
ENV PATH /usr/lib/rstudio-server/bin/:$PATH

## Download and install RStudio server & dependencies
## Attempts to get detect latest version, otherwise falls back to version given in $VER
# Step 10
RUN apt-get update \
&& apt-get install -y \
libapparmor1 \
libcurl4-openssl-dev \
psmisc \
python-setuptools

# Step 11
RUN VER=$(wget --no-check-certificate -qO- https://s3.amazonaws.com/rstudio-server/current.ver) \
&& wget -q http://download2.rstudio.org/rstudio-server-${VER}-amd64.deb \
&& dpkg -i rstudio-server-${VER}-amd64.deb \
&& rm rstudio-server-*-amd64.deb \
&& ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc /usr/local/bin \
&& ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc /usr/local/bin \
&& wget https://github.com/jgm/pandoc-templates/archive/1.15.0.6.tar.gz \
&& mkdir -p /opt/pandoc/templates && tar zxf 1.15.0.6.tar.gz \
&& cp -r pandoc-templates*/* /opt/pandoc/templates && rm -rf pandoc-templates* \
&& mkdir /root/.pandoc && ln -s /opt/pandoc/templates /root/.pandoc/templates

## Ensure that if both httr and httpuv are installed downstream, oauth 2.0 flows still work correctly.
# Step 12
RUN echo '\n\
\n# Configure httr to perform out-of-band authentication if HTTR_LOCALHOST \
\n# is not set since a redirect to localhost may not work depending upon \
\n# where this Docker container is running. \
\nif(is.na(Sys.getenv("HTTR_LOCALHOST", unset=NA))) { \
\n  options(httr_oob_default = TRUE) \
\n}' >> /etc/R/Rprofile.site

# Step 13
RUN wget -P /tmp/ https://github.com/just-containers/s6-overlay/releases/download/v1.11.0.1/s6-overlay-amd64.tar.gz \
  && tar xzf /tmp/s6-overlay-amd64.tar.gz -C /

## Clean up installations
# Step 14
RUN apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Step 15
COPY userconf.sh /etc/cont-init.d/conf
# Step 16
COPY run_rstudio_server.sh /etc/services.d/rstudio/run
# Step 17
COPY rocker_Renviron /home/rstudio/.Renviron
# Step 18
COPY rocker_Rprofile /home/rstudio/.Rprofile

# Step 19
## Install additional required R packages
## Install all required R packages
RUN  /bin/bash -c "export TM_HOME='/tissue_miner'; /tissue_miner/Setup.R | tee ${TM_HOME}/.tm_install_rsetup.log" \
	/bin/bash -c "export TM_HOME='/tissue_miner'; /tissue_miner/misc/rstudio_TM_docker_image/RstudioSetup.R | tee ${TM_HOME}/.tm_install_rstudio_setup.log"

# Step 20
# RUN ln -s /tissue_miner /home/rstudio/tissue_miner
RUN chmod -R ugo+rwx /root

# Step 21
EXPOSE 8787


## also add .bash_profile to user home directory
COPY docker_bash_profile.sh /.bash_profile

## Source .bash_profile when running the container
COPY docker_entrypoint.sh /
ENTRYPOINT ["/docker_entrypoint.sh"]

# Step 22
CMD ["/init"]


