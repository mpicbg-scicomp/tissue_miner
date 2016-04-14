

![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/stripes_0.jpg)
![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/veins_0.jpg)


About
=================

TissueMiner is a framework to perform a multiscale analysis of an entire developing epithelium acquired at cellular resolution over several hours or days. It  includes **(a)** tools to build a database from a timelapse movie, **(b)** ready-to-use tools to analyze various aspects of cell and tissue morphogenesis, and **(c)** a powerful but yet easy to use R programming interface to allow for the implementation of custom analyses.

* [Installation](https://github.com/mpicbg-scicomp/tissue_miner#installation)
* [Documentation](https://github.com/mpicbg-scicomp/tissue_miner#documentation)

### (a) Cell Tracking Data Database


First, TissueMiner tracks cells in timelapse image series of a tissue and stores the results in a database. To do so, it uses [TissueAnalyzer](MovieProcessing.md#TissueAnalyzer) to segment and track cells in time and space.

Next, TissueMiner extracts and organizes information about cell geometry, topology (cell neighbor relationships) and ancestry from the TissueAnalyzer outputs into a database. By using sqlite, each movie is essentially converted into a _single_ database file.


### (b) Tools to Analyse Epithelium Morphogenesis

TissueMiner comes along with a great variety of tools to work with the created movie database. These include
* Lineage Browsing
* Cell Shear Contribution Analysis
* Topology Analyses (cell rearrangements, cell packing, etc.)
* and some more. See the [Movie Processing Guide](MovieProcessing.md#Tools) for a complete listing

Many of these tools provide tabular results as well as rendered movies to make the complex nature of timelapse more accessible.


### (c) High-Level Programming API

TissueMiner provides a convenient R programming interface to query, quantify and visualize cell dynamics during epithelium morphogenesis in [R](https://www.r-project.org/). The first step of a custom workflow is to query the movie database to the extract properties like:

* Cell state properties (position, area, anisotropy, cell packing geometry, fluorescence intensity)
* Rates of cellular events (divisions, cell neighbor changes, extrusions, shape changes)
* Orientation of cellular events (unit nematic description)
* Rates of deformations of each type of cellular event (tensorial description)
* Rates of tissue deformation (area expansion and convergence-extension) contributed by each type of cellular event (tensorial description)

These data can than be employed to

* Visualize of quantified data in graphs or directly on the movie images
* Do statistics (time evolution of the distributions of cell area, anisotropy, packing, ..., bond length, vertex multiplicity, ...)
* Synchronize different movies in time
* Comparing values between multiple movies and ROI's
* Do multiscale quantification and visualization using both tracked ROI's (Lagrangian description) and fixed grids (Eulerian description)



Installation
================

TissueMiner ships with a command-line Ubuntu-Installer and can be used via Docker on other platforms.

### Ubuntu

* To install it on Ubuntu simply clone this repository and run the installation script `install_tm.sh` in a [terminal](https://help.ubuntu.com/community/UsingTheTerminal) as follow:

```
 ## Define the directory where to install TissueMiner
 export TM_HOME="${HOME}/tissue_miner"

 ## Download TissueMiner repository
 sudo apt-get install git
 git clone https://github.com/mpicbg-scicomp/tissue_miner.git ${TM_HOME}
 
 ## Install TM
 ${TM_HOME}/installation/ubuntu/install_tm.sh
 source ~/.bashrc
```

* To update it on Ubuntu, click [here](faq.md#how-to-update-my-tissueminer-installation)

### Other Systems

* To install TissueMiner on MacOS, Windows or any non-Ubuntu Linux system we provide a Docker container that bundles TissueMiner and all its dependencies. If not yet present on your system, you need to install the [docker engine](https://docs.docker.com/)
beforehand.

* Next, you can download the TissueMiner application bundled in a docker image called _etournay/tissue_miner_.

* On Mac or Windows, open a Docker Quick Start Terminal: ![alt tag](docs/readme_screenshots/docker_toolbox_osx_quickstart_icon_nolabel.png)
```
docker pull etournay/tissue_miner
```

* For troubleshooting in case of out-of-date operating systems [click here](misc/docker_troubleshooting.md).

* To update it, click [here](faq.md#how-to-update-my-tissueminer-installation)

Documentation
================

* To get started with TissueMiner we suggest to do **Quickstart Tutorials:**
    + **[Quickstart R Tutorials](docs/quickstart/ubuntu/tm_qs_example_data.md#first-use-of-tissueminer-with-example-data)** (Ubuntu).
    + **[Quickstart docker-R Tutorials](docs/quickstart/other_os/tm_qs_example_data.md#first-use-of-tissueminer-with-example-data)** (any systems).

* To get started with your own data analysis:
    + **[R Tutorials](docs/quickstart/ubuntu/tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)** (Ubuntu).
    + **[docker-R Tutorials](docs/quickstart/other_os/tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)** (any systems).

* We also provide an exhaustive **[TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html)** with lots of examples, background information and API details:
    + learn the necessary basics of the R language for TissueMiner
    + learn the powerful libraries of TissueMiner for epithelium analysis
    + benefit from lots of examples to make your own scripts
    
* You already know Python and you like to use it instead ? We also provide a **[Tutorial in Python](docs/TM_tutorial_in_Python/TissueMiner_pythonTutorial-3WT_Demo.md#tissueminer-python-tutorial)**

* For a more general overview consider our [Resource Paper](/link/here/once/published).

* Commonly asked questions are answered in the [FAQ section](faq.md).

Datasets
================
We provide four datasets:
* 1 demo sample (~150 Mb) to run [Quickstart R Tutorials](docs/quickstart/tm_quickstart_landing_page.md#first-use-of-tissueminer-from-the-command-line)
* 3 big datasets (>2Gb each) to explore more advanced capabilities of TissueMiner described in the [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html)

Support
=========

Please use the [github ticket system](https://github.com/mpicbg-scicomp/tissue_miner/issues) to report issues or suggestions. We also welcome pull requests.



Reference
==========

If you like to use TissuMiner for your own research, please cite

> Etournay et al. (2015). Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife, 4, e07090. [doi:10.7554/eLife.07090](http://elifesciences.org/content/early/2015/06/23/eLife.07090)


