

![alt tag](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/readme_screenshots/stripes_0.jpg)
![alt tag](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/readme_screenshots/veins_0.jpg)

News
=================

**New TM release on Nov 18th:**

* We now highly recommend to run TM through Docker (Win7, Win10, MacOS, other Unix)
* Designed to work through Docker to simplify installation and maintenance
* Versionning now enabled for docker images: new release is 1.0.1
* Please update the tm alias in your .bashrc:

```
alias tm='docker run --rm -ti -v $(dirname $PWD):/movies -w /movies/$(basename $PWD) -u rstudio etournay/tissue_miner'
```



About
=================

TissueMiner is a framework to perform a multiscale analysis of an entire developing epithelium acquired at cellular resolution over several hours or days. It  includes **(a)** tools to build a database from a timelapse movie, **(b)** ready-to-use tools to analyze various aspects of cell and tissue morphogenesis, and **(c)** a powerful but yet easy to use R programming interface to allow for the implementation of custom analyses.

* [Installation](https://github.com/mpicbg-scicomp/tissue_miner#installation)
* [Documentation](https://github.com/mpicbg-scicomp/tissue_miner#documentation)

### (a) Cell Tracking Data Database

First, TissueMiner tracks cells in timelapse image series of a tissue and stores the results in a database. To do so, it uses [TissueAnalyzer](MovieProcessing.md#TissueAnalyzer) to segment and track cells in time and space.

Next, TissueMiner extracts and organizes information about cell geometry, topology (cell neighbor relationships) and ancestry from the TissueAnalyzer outputs into a database. By using sqlite, each movie is converted into a _single_ database file.


### (b) Tools to Analyse Epithelium Morphogenesis

TissueMiner comes along with a great variety of tools to work with the created movie database. These tools include
* Cell geometry analysis
* Cell topology analysis
* Cell lineage browsing
* Cell contributions to tissue deformation

Many of these tools provide tabular results as well as rendered movies to make the complex nature of timelapse more accessible.

### (c) High-Level Application Programming Interface (API)

TissueMiner provides a convenient R programming interface to query, quantify and visualize cell dynamics during epithelium morphogenesis in [R](https://www.r-project.org/). The functions provided with the TissueMiner API allow one to query the movie database for extracting quantities like:

* Cell state properties (position, area, anisotropy, cell packing geometry, fluorescence intensity)
* Rates and orientation of cellular events (divisions, cell neighbor changes, extrusions, shape changes)
* Rates of tissue deformation (area expansion and convergence-extension) contributed by each type of cellular event (tensorial description)

This API can then be employed to

* Visualize quantified data in graphs or directly on the movie images
* Do statistics (distribution of cell area, anisotropy, packing, bond length, vertex multiplicity, ...)
* Synchronize different movies in time
* Compare values between multiple movies and ROI's
* Do multiscale quantification and visualization using both tracked ROI's (Lagrangian description) and fixed grids (Eulerian description)



Installation
================

TissueMiner ships with a command-line Ubuntu-Installer and can be used via Docker on other platforms.

### Ubuntu

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

* To update TissueMiner on Ubuntu, click [here](faq.md#how-to-update-my-tissueminer-installation)

### Other Systems

* To install TissueMiner on a MacOS or a Windows system, we provide a Docker container that bundles TissueMiner and all its dependencies. If not yet present on your system, you need to install the [docker toolbox](https://www.docker.com/products/docker-toolbox)
beforehand (180Mb). On any non-Ubuntu Linux, please install the [docker engine](https://docs.docker.com/).

* Next, you can download the TissueMiner application bundled in a docker image called _etournay/tissue_miner_ (~700Mb).

    + On Mac or Windows, open a Docker Quick Start Terminal: ![alt tag](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/readme_screenshots/docker_toolbox_osx_quickstart_icon_nolabel.png)
```
docker pull etournay/tissue_miner
```

* For troubleshooting in case of out-of-date operating systems [click here](misc/docker_troubleshooting.md).

* To update TissueMiner, click [here](faq.md#how-to-update-my-tissueminer-installation)

Documentation
================

* First [install](https://github.com/mpicbg-scicomp/tissue_miner#installation) TissueMiner on your system.
* To get started with TissueMiner we suggest to do **Quickstart Tutorials on our example data:** 
    + **[Quickstart TM-R Tutorials](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/quickstart_tutorial/ubuntu/tm_qs_example_data.md#first-use-of-tissueminer-with-example-data)** (Ubuntu).
    + **[Quickstart docker-TM-R Tutorials](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/quickstart_tutorial/other_os/tm_qs_example_data.md#first-use-of-tissueminer-with-example-data)** (Mac, Windows, other Linux).

* To get started with your **own data** analysis:
    + **[TM-R Tutorials](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/quickstart_tutorial/ubuntu/tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)** (Ubuntu).
    + **[docker-TM-R Tutorials](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/quickstart_tutorial/other_os/tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)** (Mac, Windows, other Linux).
    
* We provide the streamlined [R scripts](docs/quickstart/scripts) that are used in the tutorials to perform a simple analysis of a single movie. This includes a script entitled `analyse_movie.R` that performs most of the analyses described in the Quickstart Tutorials.

* We also provide an exhaustive **[TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html)** with lots of examples, background information and API details:
    + learn the necessary basics of the R language for TissueMiner (also see [here](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/Learning_the_R_basics_for_TissueMiner.html))
    + learn the powerful libraries of TissueMiner for epithelium analysis
    + benefit from lots of examples to make your own scripts
    
* You already know Python and you like to use it instead ? We also provide a **[Tutorial in Python](docs/TM_tutorial_in_Python/TissueMiner_pythonTutorial-3WT_Demo.md#tissueminer-python-tutorial)**

* For a more general overview consider our [Resource Paper](http://dx.doi.org/10.7554/eLife.14334).

* Commonly asked questions are answered in the [FAQ section](faq.md).

Datasets
================
We provide four datasets:
* 1 demo sample (~35 Mb) to run [Quickstart TM-R Tutorials](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/quickstart_tutorial/tm_quickstart_landing_page.md#first-use-of-tissueminer)

```{bash}
# In your unix terminal (or Docker QuickStart Terminal), type in:
curl https://cloud.mpi-cbg.de/index.php/s/vC0VqD2Wy4A6uqu/download  | tar -zxvf -
```

* 3 big datasets (~800Mb each) to explore more advanced capabilities of TissueMiner described in the [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html)

```{bash}
# In your unix terminal (or Docker QuickStart Terminal), type in:

# Dataset WT_1 (~800Mb)
curl https://cloud.mpi-cbg.de/index.php/s/SgxxQk5CkIpTLPW/download  | tar -zxvf -

# Dataset WT_2 (~800Mb)
curl https://cloud.mpi-cbg.de/index.php/s/Z6ZR1b0sGWnC8Cj/download  | tar -zxvf -

# Dataset WT_3 (~600Mb)
curl https://cloud.mpi-cbg.de/index.php/s/4BJiyxnCS7HFyKB/download  | tar -zxvf -
```


Support
=========

Please use the [github ticket system](https://github.com/mpicbg-scicomp/tissue_miner/issues) to report issues or suggestions. We also welcome pull requests.



Reference
==========

If you like to use TissuMiner for your own research, please cite

> Etournay et al. (2016). TissueMiner: a multiscale analysis toolkit to quantify how cellular processes create tissue dynamics. eLife, 5, e14334. [doi:10.7554/eLife.14334](https://elifesciences.org/content/5/e14334)

or

> Etournay et al. (2015). Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife, 4, e07090. [doi:10.7554/eLife.07090](http://elifesciences.org/content/early/2015/06/23/eLife.07090)



