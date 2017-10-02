#### CHANGES IN VERSION 1.x.x

* Release date: **TBA**
* New shear implementation using infinitesimal displacement formalism (Merkel et al 2016)
* New mqf* libraries that enable interpolation of data, keep the old library still available
* New Docker image including newer version of R, RstudioServer and dependencies
* Host user ID mapping to the docker container

#### CHANGES IN VERSION 1.0.3

* Release date: **October 9th, 2017**
* Bug fix in image format conversion by imagemagick (creation of original.png from RGB files)
* Bug fix in applying user-defined configuration files

#### CHANGES IN VERSION 1.0.2

* Release date: **May 1st, 2017**
* Bug fixe in UserRoiTracking.R: fixed border cell status ambiguities that caused duplicated cell_id and merge failure in rare cases
* Commited changes in the existing docker image

#### CHANGES IN VERSION 1.0.1

* Release date: **Nov 16th, 2016**
* Enabled version check for updates
* Security: now runs the workflow with user permissions and solve permission issues when using docker on Linux
* Control of resource allocation
* Can now run the workflow without network connection
* Corresponding [user manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual_v1.0.1.html)

#### CHANGES IN VERSION 1.0.0

* Release date: **May 2016**
* RStudio server was added to the docker image for running TissueMiner API from within a fully configured RStudio and R environment along with a fully configured Ubuntu environment
* Corresponding [user manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html)

