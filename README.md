

![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/stripes_0.jpg)
![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/veins_0.jpg)


About
=================

TissueMiner is a framework to perform a multiscale analysis of an entire developing epithelium acquired at cellular resolution over several hours or days. It  includes tools to **(a)** build a database from timelapse movie, **(b)** ready-to-use tools to analyze various aspects of cell morphogenesis, and **(c)** a powerful but yet easy to use R programming interface to allow for the implementation of custom analyses.

### Cell Tracking Data Database


First, TissueMiner tracks cells in timelapse image series of tissues and stores the results in a database. To do so, it uses [TissueAnalyzer](MovieProcessing.md#TissueAnalyzer) to segment and track cells in time and space.

Next, TissueMiner extracts and organizes information about cell geometry, topology (cell neighbor relationships) and ancestry from the TissueAnalyzer outputs into a database. By using sqlite each movie is essentially converted into a _single_ database file.


### Tools to Analyse Epithelium Morphogenesis

TissueMiner comes along with a great variety of tools to work with the created movie database. These include
* Lineage Tracking
* Cell Shear Contribution Analysis
* Topolgy Anlayses (T1s, polygon classes, etc.)
* and some more. See the [Batch Processing Guide](MovieProcessing.md#Tools) for a complete listing

Many of these tools provide tabular results as well as rendered movies to make the complex nature of timelapse more accessible.


TissueMiner also contains a GUI (**todo rapha** really where is it?) to group cells (ROI's) in space, and a lineage browsing algorithm to follow these ROI's backward and forward in time.


### High-Level Programming API

**todo rapha** mayebe we could thin out this section and put the complete listings into the R manual?

TissueMiner provides a convenient interface to query, quantify and visualize cell dynamics during epithelium morphogenesis in [R](https://www.r-project.org/).

By doing so the following properties can be accessed
* Cell state properties (position, area, anisotropy, cell packing geometry, fluorescence intensity)
* Rates of cellular events (divisions, cell neighbor changes, extrusions, shape changes)
* Orientation of cellular events (unit nematic description)
* Rates of deformations of each type of cellular event (tensorial description)
* Rates of tissue deformation (area expansion and convergence-extension) contributed by each type of cellular event (tensorial description)

It provides means to
* Analyze time and tissue-orientation registrations for comparing multiple movies
* Synchronize different movies
* Compare tissue and cell deformations between ROI's and between movies.
* do statistics (time evolution of the distributions of cell area, anisotropy, packing, ..., bond length, vertex multiplicity, ...)
* Do multiscale quantification and visualization using both dynamic ROI's (Lagrangian description) and fixed grids (Eulerian description): from individual cells to averages over the entire tissue
* Visualize of quantified data in graphs or directly on the movie images


Installation
================

TissueMiner ships with an Ubuntu-Installer and can be used via Docker on other platforms.

### Ubuntu

To install it on Ubuntu simply clone this repository and run the installation script `install_tm.sh`

```
    ## Define the directory where to install TissueMiner
    export TM_HOME="~/tissue_miner"

    ## download this repository
    git clone https://github.com/mpicbg-scicomp/tissue_miner.git ${TM_HOME}

    sudo ${TM_HOME}/install_tm.sh | tee ${TM_HOME}/installation.log
```

To test the installation run
```
${TM_HOME}/tm_test.R
```

### Other Systems

To install TissueMiner on MacOS, Windows or any non-Ubuntu Linux system we provide a Docker container that bundles TissueMiner and all its dependencies. If not yet present on your system, you need to install the [docker engine](https://docs.docker.com/)
beforehand.


Next, you can download the TissueMiner application bundled in a docker image called _brandl/tissue_miner_
```
docker pull brandl/tissue_miner
```

To test the installation with the following command:
```
docker run --rm -ti -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner /bin/bash ${TM_HOME}/tm_test.R
```

For troubleshooting in case of out-of-date operating systems [click here](misc/docker_troubleshooting.md).



Documentation
================

To get started with TissueMiner we suggest to do the **[quickstart tutorial](https://mpicbg-scicomp.github.io/tissue_miner/quickstart/TM_Quickstart.html)** using provided example data.

We also provide an exhaustive [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html) with lots of examples, background information and API details.

To learn about how to process your own movies see the [Batch Processing Guide](MovieProcessing.md).

For a more general overview consider our [Resource Paper](/link/here/once/published).

Commonly asked questions are answered in the [FAQ section](faq.md).


Support
=========

Please use the [github ticket system](https://github.com/mpicbg-scicomp/tissue_miner/issues) to report issues or suggestions. We also welcome pull requests. If



Reference
==========

If you like to use TissuMiner for your own research, please cite

> Etournay et al. (2015). Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife, 4, e07090. [doi:10.7554/eLife.07090](http://elifesciences.org/content/early/2015/06/23/eLife.07090)


