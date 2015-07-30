#######################################################
### Build Image using Dockerfile
#######################################################

export TM_HOME=/home/brandl/mnt/mack/project-raphael/scripts/tissue_miner
cd ${TM_HOME}/misc

docker build -t brandl/tissue_miner .

## publish the new image
#docker login
#docker push brandl/tissue_miner

docker run -t -i -v $(pwd)/example_movies:/movies brandl/tissue_miner /bin/bash --login -c "cd demo_ForTissueMiner; sm -j 5 all"

docker run -t -i brandl/tissue_miner /bin/bash --login


#######################################################
### Start interactive docker session (skip for local installation)
#######################################################

# https://docs.docker.com/articles/basics/
docker pull ubuntu

docker run -i -t ubuntu /bin/bash

## run all following bits within docker shell

# source setup.sh


#######################################################
## Install all docker specific stuff
#######################################################


## prepare xvfb to allow to run it without X
# run without x # run without xhttps://linuxmeerkat.wordpress.com/2014/10/17/running-a-gui-application-in-a-docker-container/
apt-get install --assume-yes xvfb

#Xvfb :1 -screen 0 1024x768x16 &> xvfb.log  &
#DISPLAY=:1.0

#######################################################
## Commit to local docker registry (skip for local installation)
#######################################################

## submit changes to local docker registry
## http://stackoverflow.com/questions/19585028/i-lose-my-data-when-the-container-exits


#sudo docker ps -l
#sudo docker commit a8889839969b brandl/tissue_miner
#sudo docker commit a38860f9a922 brandl/test
sudo docker commit $(sudo docker ps -l | cut -f1 -d' ' | tail -n+2) brandl/tissue_miner

#sudo docker run iman/ping ping www.google.com
docker run -i -t brandl/tissue_miner /bin/bash
#docker run -i -t brandl/test /bin/bash

## export docker image for testing on other machine
#http://www.jamescoyle.net/how-to/1512-export-and-import-a-docker-image-between-nodes
docker save brandl/tissue_miner > ~/mnt/mack/project-raphael/tissue_miner.docker.$(date +'%y%m%d').tar

## dockerfile
#http://stackoverflow.com/questions/19585028/i-lose-my-data-when-the-container-exits
# --> Incrementally committing changes is not "the docker way". Use a DOCKERFILE. –  user2105103 Jan 14 at 14:23
# https://www.digitalocean.com/community/tutorials/docker-explained-using-dockerfiles-to-automate-building-of-images


#docker ps -a
#docker start 6f0d6b6a503b
#docker attach 5fa0912c33ee



#######################################################
## Run the workflow via docker
#######################################################

## grab a slim copy
#cd ~/projects/flywing
#cp -r ~/mnt/mack/project-raphael/movie_dbs/db_tests/tissue_analyzer/demo_ForTissueMiner ~/projects/flywing/

## run sm wrapper interactivly and as tool
#docker run -t -i -v /bin/bash --login
docker run -t -i brandl/tissue_miner /bin/bash --login -c "source .bash_profile; sm --help"
docker run -t -i brandl/tissue_miner /bin/bash --login -c "source .bash_profile; sm --help"
docker run -t -i brandl/tissue_miner /bin/bash --login -c "sm --help"
docker run -t -i brandl/tissue_miner /bin/bash --login


docker run -t -i -v /Users/brandl/Desktop/my_movies:/movies brandl/tissue_miner /bin/bash --login -c "source .bash_profile; cd demo_ForTissueMiner; sm -n all"
docker run -t -i -v /Users/brandl/Desktop/my_movies:/movies brandl/tissue_miner /bin/bash --login -c "source .bash_profile; cd demo_ForTissueMiner; sm makedb"
docker run -t -i -v /Users/brandl/Desktop/my_movies:/movies brandl/tissue_miner /bin/bash --login
#docker run -t -i -v /home/brandl/projects/flywing/:/movies brandl/tissue_miner /bin/bash --login -c "cd demo_ForTissueMiner; sm  --force -p area_movies"
#docker run -t -i -v /home/brandl/projects/flywing/:/movies brandl/tissue_miner /bin/bash --login -c "cd demo_ForTissueMiner; sm --force makedb"

## docker cheatsheet https://github.com/wsargent/docker-cheat-sheet


#######################################################
## Test TissueMiner within the interactive docker shell
#######################################################


movieName="$(basename $(pwd))"

sm -n -p

sm -p --force roi_tracking
sm -p --force shear_calculate

## try the parsering
sm -p parse_tables

sm -p all | tee ${movieName}.sm.log

## print some statistics
sm -S > sm_summary.txt
sm  --dag | dot -Tpdf > dag_tbd.pdf

#mailme "${movieName}: rapha done"

## trash results and start over again (double disabled command to prevent accidental usage)
#find . -type f -not -name 'Segmentation' #| xargs rm -rf

### segfault debugging

# downgrade dplyr http://stackoverflow.com/questions/31010713/combining-dplyrmutate-with-lubridateymd-hms-in-r-randomly-causes-segfault


#######################################################
## Misc
#######################################################

http://askubuntu.com/questions/472412/how-do-i-upgrade-docker

## bundle the example movies
tar -zcvf tm_example_data.tar.gz example_movies