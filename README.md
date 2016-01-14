
About
=================

**TissueMiner** is a framework to track cells and groups of cells (rois) in space an and time. It allows to calculate statistics about lineage, shape, deformation and topology. It comes along with various tools to visualize these data. 


![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/stripes_0.jpg)
![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/veins_0.jpg)


How to get started?
================
Install the [docker engine](https://docs.docker.com/) on your system.

We provide a preconfigured [docker image](https://registry.hub.docker.com/u/brandl/tissue_miner/) to run TissueMiner without any setup. Just do

     ## download some example data
     curl https://cloud.mpi-cbg.de/index.php/s/EspCWSQn3k6NKLA/download  | tar -zxvf -
     
     ## download the image and start the analysis
     docker pull brandl/tissue_miner
     docker run --rm -ti -v $(pwd)/example_movies:/movies -w /movies/demo_ForTissueMiner brandl/tissue_miner sm all
     
To run TissueMiner over your own data, you'll need to set the source directory containing your movie directories and to replace the example movie name with your movie/directory of interest.


How to run locally?
================

First, make sure that the packages listed in [install_dependencies.sh](misc/install_dependencies.sh) are installed on your system.

Second, just grab a copy of TissueMiner and run the setup procedure in a terminal (Unix-shell interpreter)

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

Then you can run TissueMiner to analyze your movie, provided your movie folder contains a "Segmentation" folder in which to find the original images as well as the TissueAnalyzer outputs.

To run the workflow we recommend [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home). We provide a [snakemake workflow](workflow/tm.snkmk) to ease running TissueMiner on a cluster or locally on a single computer. It integrates all analyses implemented in TissueMiner and can be easily extended to include project specifc elements as well. This is how we usually run TissueMiner: 

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

Data structure
================

We advice the user to store all movies in a movie repository folder <movie_repository> to facilitate automated movie comparison.

Here is the required structure of a movie:


\<movie_repository\>/\<movie_directory\>/**Segmentation**/\<movie_directory_name\>_%03d.png

%03d represents a number padded with 3 digits. The number of digits can be modified [FAQ](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/faq.md).

All movie images are contained in the 'Segmentation' folder.
Upon segmentation and tracking, TissueAnalyzer will generate one folder per image, in the Segmentation folder.

TissueMiner will generate additional folders and files in the \<movie_directory\> folder.

Documentation
=============

We provide an exhaustive **[tutorial](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html)** with lots of examples, background information and API details. For a more general overview consider the [resource paper](/link/here/once/published)

TissueAnalyzer
=============
-Tissue Analyzer (copyright Aigouy 2016) ships with its own licence
(see license_TA.txt bundled in the software). Tissue Analyzer should not be modified or
reverse engineered.  Tissue Analyzer should always be distributed
bundled with TM and not alone.

You can install the latest version of Tissue Analyzer (formerly known
as Packing Analyzer) within FIJI (http://fiji.sc/Fiji). To do so:

* -launch FIJI
* -open the "Help" menu
* -click on "Update..."
* -click on "Manage update sites"
* -click on "Add" and enter "http://sites.imagej.net/TA/" in the URL field (and anything you like in the "Name" field)
* -click "Close" and FIJI should offer you to install TA.

Once the installation is complete, restart FIJI, open the "Plugins" menu and click on "Tissue Analyzer"

For a quick start guide, click [here](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/docs/TAdoc.pdf)

Support
=========

Please use the github ticket system to report issues or suggestions. We also welcome pull requests.


Reference
==========

If you like to use TissuMiner for your own research, please cite

> Etournay et al. (2015). Interplay of cell dynamics and epithelial tension during morphogenesis of the Drosophila pupal wing. eLife, 4, e07090. [doi:10.7554/eLife.07090](http://elifesciences.org/content/early/2015/06/23/eLife.07090)


