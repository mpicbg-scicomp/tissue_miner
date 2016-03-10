

![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/stripes_0.jpg)
![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/veins_0.jpg)


About
=================

**TissueMiner** is a framework to perform a multiscale analysis of an entire developing epithelium acquired at cellular resolution over several hours or days. It works in [three steps](misc/TM_description.md). 


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
on your system.


Next, we can download the TissueMiner application bundled in a docker image called _brandl/tissue_miner_
```
docker pull brandl/tissue_miner
```

To test the docker installation you can print the tissue miner version:
```
docker run --rm -ti -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner /bin/bash ${TM_HOME}/tm_test.R
```

For troubleshooting in case of out-of-date operating systems [click here](misc/docker_troubleshooting.md).



Documentation
================

To get started with TissueMiner we suggest to do the **[quickstart tutorial](https://mpicbg-scicomp.github.io/tissue_miner/quickstart/TM_Quickstart.html)** using provided example data.

We aslo provide an exhaustive [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html)** with lots of examples, background information and API details. For a more general overview consider the [resource paper](/link/here/once/published)

To learn about how to process your own movies see the [Batch Processing Guide](MovieProcessing.md)


Support
=========

Please use the [github ticket system](https://github.com/mpicbg-scicomp/tissue_miner/issues) to report issues or suggestions. We also welcome pull requests.


FAQ
=========
Go to [the FAQ page](faq.md)

Reference
==========

If you like to use TissuMiner for your own research, please cite

> Etournay et al. (2015). Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife, 4, e07090. [doi:10.7554/eLife.07090](http://elifesciences.org/content/early/2015/06/23/eLife.07090)


