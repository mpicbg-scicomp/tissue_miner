

![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/stripes_0.jpg)
![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/veins_0.jpg)


About
=================

**TissueMiner** is a framework to perform a multiscale analysis of an entire developing epithelium acquired at cellular resolution over several hours or days. It works in [three steps](misc/TM_description.md). 


How to get started?
================
Install the [docker engine](https://docs.docker.com/) on your system. And start it (a shell window pops up). Docker can run the TissueMiner automated workflow on your system (MacOSX, Linux, Windows) without any OS setup. To do so, you once need to download the TissuerMiner [docker image](https://registry.hub.docker.com/u/brandl/tissue_miner/) that we preconfigured. Just use the opened shell window (Unix-shell interpreter) and do

     ## DEMO:
     ## Step 1: download some example data for the demo
     curl https://cloud.mpi-cbg.de/index.php/s/EspCWSQn3k6NKLA/download  | tar -zxvf -
     
     ## Step 2: download the TissueMiner application (bundled in a docker image called "brandl/tissue_miner")
     docker pull brandl/tissue_miner
     
     ## Step 3: Run the TissueMiner automated workflow on the provided example
     docker run --rm -ti -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner sm all

For troubleshooting in case of out-of-date operating systems [click here](misc/docker_troubleshooting.md).

To run TissueMiner over **your own data**, you'll need to go to the source directory containing your movie directories (use the *cd* command + drag and drop the folder containing your movie directories onto the shell + press enter) and to replace the example movie name (demo_ForTissueMiner) with your movie/directory of interest.

    ## Example for "my_favorite_movie"
    cd path_to_your_movie_repository
    docker run --rm -ti -v $(pwd):/movies -w /movies/my_favorite_movie brandl/tissue_miner sm all

Once your tracked-cell data have been processed by the TissueMiner automated workflow, you can perform a custom multiscale analysis of epithelial morphogenesis using the detailed [R-tutorial](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html) or the [Python-tutorial](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/docs/TM_tutorial_in_Python/TissueMiner_pythonTutorial-3WT_Demo.md).


Data structure
================

We advice the user to store all movies in a movie repository folder \<movie_repository\> to facilitate automated movie comparison.

Here is the **required structure** of a movie:


\<movie_repository\>/\<movie_directory\>/**Segmentation**/\<movie_directory_name\>_%03d.png

%03d represents a number padded with 3 digits. The number of digits can be modified [FAQ](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/faq.md).

All movie images are contained in the **Segmentation** folder.
Upon segmentation and tracking, TissueAnalyzer will generate one folder per image, in the Segmentation folder.

TissueMiner will generate additional folders and files in the \<movie_directory\> folder.

Installation on Ubuntu
================

First, make sure that the packages listed in [install_dependencies.sh](misc/install_dependencies.sh) are installed on your Ubuntu system.

Second, run the setup procedure in a terminal (Unix-shell interpreter)

    ## Set path to install TissueMiner in you home folder (please use an absolute path)
    export TM_HOME="~/tissue_miner"

    ## download this repository
    git clone https://github.com/mpicbg-scicomp/tissue_miner.git ${TM_HOME}
        
    ## or update existing local copy with
    cd  ${TM_HOME} && git pull origin && cd parser && make clean all
    
    ## Install all required R packages (or run Setup.R in RStudio)
    ${TM_HOME}/Setup.R | tee tm_setup.log
    
    ## compile the parser needed to convert TissueAnalyzer outputs into csv
    cd ${TM_HOME}/parser && make

    ## adjust your path to include all tools
    export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
    export PATH=${TM_HOME}/parser:$PATH
    

Don't forget to define a TM_HOME shell variable, pointing to the root of your TissueMiner installation, since it will require it to resolve script paths internally.

Then you can run TissueMiner to analyze your movie, provided your movie folder contains a **Segmentation** folder in which to find the original images as well as the TissueAnalyzer outputs.

To run the automated workflow we recommend [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home). We provide a [snakemake workflow](workflow/tm.snkmk) to ease running TissueMiner on a cluster or locally on a single computer. It integrates all analyses implemented in TissueMiner and can be easily extended to include project specifc elements as well. This is how we usually run TissueMiner: 

    ##  define a custom snakemake launcher to save typing
    sm() {
        snakemake --snakefile ${TM_HOME}/workflow/tm.snkmk --keep-going "$@"
    }
    export -f sm
    
    ## TissueMiner assumes the movie data to be present in the current working directoy
    cd <movie_directory>
    
    ## List all tasks
    sm -n
    
    ## Process all tasks and write a log file
    sm all | tee log.txt
    
    ## ... or just run sepecific tasks
    sm make_db | tee log.txt
    
    ## Export statistics and execution state graph visualization
    sm --dag | dot -Tpdf > dag_tbd.pdf
    sm -D > sm_execution_state.txt
    
For snakemake details see the its [reference](https://bitbucket.org/johanneskoester/snakemake/wiki/Home).

Although we do not recommend it, you can also run each of the TissueMiner tools separately. See [simple_workflow.sh](workflow/simple_workflow.sh) for an example pipeline.


We recommend to add the following lines in your *.bashrc* or *.bash_profile* so that you don't have to export the environment variables again:

    ## Copy these lines into your .bashrc (or .bash_profile):
    export TM_HOME="~/tissue_miner"
    export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
    export PATH=${TM_HOME}/parser:$PATH
    sm() {
        snakemake --snakefile ${TM_HOME}/workflow/tm.snkmk --keep-going "$@"
    }
    export -f sm
    
    ## Then apply changes
    source .bashrc



Documentation
=============

We provide an exhaustive **[tutorial](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html)** with lots of examples, background information and API details. For a more general overview consider the [resource paper](/link/here/once/published)

TissueAnalyzer
=============
**Tissue Analyzer** (copyright Aigouy 2016) ships with its own licence
(see license_TA.txt bundled in the software). Tissue Analyzer should not be modified or
reverse engineered.  Tissue Analyzer should always be distributed
bundled with **TissueMiner** and not alone.

You can install the latest version of Tissue Analyzer (formerly known
as Packing Analyzer) within FIJI (http://fiji.sc/Fiji). To do so:

* launch FIJI
* open the "Help" menu
* click on "Update..."
* click on "Manage update sites"
* click on "Add" and enter "http://sites.imagej.net/TA/" in the URL field (and anything you like in the "Name" field)
* click "Close" and FIJI should offer you to install TA.

Once the installation is complete, restart FIJI, open the "Plugins" menu and click on "Tissue Analyzer"

For a quick start guide, click [here](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/docs/TAdoc.pdf)

Support
=========

Please use the [github ticket system](https://github.com/mpicbg-scicomp/tissue_miner/issues) to report issues or suggestions. We also welcome pull requests.


Reference
==========

If you like to use TissuMiner for your own research, please cite

> Etournay et al. (2015). Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife, 4, e07090. [doi:10.7554/eLife.07090](http://elifesciences.org/content/early/2015/06/23/eLife.07090)


