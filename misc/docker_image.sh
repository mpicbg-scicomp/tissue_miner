#######################################################
### Build Image using Dockerfile
#######################################################

export TM_HOME=/home/brandl/mnt/mack/project-raphael/scripts/tissue_miner
cd ${TM_HOME}/misc || exit

docker build -t brandl/tissue_miner .
mailme "docker build done"

## publish the new image
#docker login
#docker push brandl/tissue_miner


## test the example from the README.md
# mcdir /home/brandl/projects/tm_demo
# curl https://cloud.mpi-cbg.de/index.php/s/EspCWSQn3k6NKLA/download  | tar -zxvf -

docker run --rm -ti -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner sm all

docker run --rm -ti -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner /bin/bash --login


#docker run --rm -ti -v $(pwd)/example_movies:/movies brandl/tissue_miner /bin/bash --login
#docker run --rm -ti -v $(pwd)/example_movies:/movies brandl/tissue_miner /bin/bash --login -c "cd demo_ForTissueMiner; sm -j 5 all"
#docker run --rm -ti -v $(pwd)/example_movies:/movies brandl/tissue_miner /bin/bash --login -c "source /.bash_profile; cd demo_ForTissueMiner; sm all"
docker run --rm -t -i -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner sm -n all
docker run --rm -t -i -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner /bin/bash --login


docker run --rm -ti -v $(pwd):/movies brandl/tissue_miner /bin/bash --login
docker run -ti -v $(pwd):/movies brandl/tissue_miner /bin/bash --login
docker run -ti -v $(pwd):/example_movies brandl/tissue_miner /bin/bash
docker run -i -t  brandl/tissue_miner /bin/bash

## signal forwarding
docker run --rm -ti -v $(pwd):/movies brandl/tissue_miner  /bin/bash --login -c "sleep 10000"
docker run --rm -ti -v $(pwd):/movies brandl/tissue_miner ls
docker run --rm -ti -v $(pwd):/movies brandl/tissue_miner /bin/bash --login
docker run --rm -ti -v $(pwd):/movies brandl/tissue_miner sm -n



#cd /Volumes/projects/project-raphael/scripts/tissue_miner
#git update-index --chmod=+x db/movie_rotation/transform_images.sh db/movie_rotation/RotateOriginals.sh

#######################################################
### Start interactive docker session (skip for local installation)
#######################################################

# https://docs.docker.com/articles/basics/
docker pull ubuntu

docker run -it ubuntu /bin/bash


#######################################################
## Commit to local docker registry (skip for local installation)
#######################################################

#sudo docker ps -l

## submit changes to local docker registry
## http://stackoverflow.com/questions/19585028/i-lose-my-data-when-the-container-exits


docker run ubuntu apt-get install -y ping

#sudo docker commit a8889839969b brandl/tissue_miner
#sudo docker commit a38860f9a922 brandl/test
docker commit $(docker ps -l | cut -f1 -d' ' | tail -n+2) brandl/test

#sudo docker run iman/ping ping www.google.com
docker run -i -t brandl/test /bin/bash
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
## TM docker debugging
#######################################################

docker run -it ubuntu /bin/bash


## do linux dependencies

docker commit $(docker ps -l | cut -f1 -d' ' | tail -n+2) brandl/test
docker run -i -t brandl/test /bin/bash

## do R packages


docker commit $(docker ps -l | cut -f1 -d' ' | tail -n+2) brandl/test
docker run -i -t brandl/test /bin/bash
## finalize image


docker commit $(docker ps -l | cut -f1 -d' ' | tail -n+2) brandl/test
docker run -i -t brandl/test /bin/bash

## test tutorial and examples
mcdir ~/projects/tm_test
curl https://cloud.mpi-cbg.de/index.php/s/EspCWSQn3k6NKLA/download  | tar -zxvf -

docker run -i -t brandl/test /bin/bash

## todo

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
docker run -t -i -v /Users/brandl/Desktop/my_movies:/movies brandl/tissue_miner /bin/bash --login -c "source .bash_profile; cd demo_ForTissueMiner; sm make_db"
docker run -t -i -v /Users/brandl/Desktop/my_movies:/movies brandl/tissue_miner /bin/bash --login
#docker run -t -i -v /home/brandl/projects/flywing/:/movies brandl/tissue_miner /bin/bash --login -c "cd demo_ForTissueMiner; sm  --force -p area_movies"
#docker run -t -i -v /home/brandl/projects/flywing/:/movies brandl/tissue_miner /bin/bash --login -c "cd demo_ForTissueMiner; sm --force make_db"

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
cd /Users/brandl/owncloud/public/tissue_miner
tar -zcvf tm_example_data.tar.gz example_movies