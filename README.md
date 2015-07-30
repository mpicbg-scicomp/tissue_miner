
About
=================

**TissueMiner** is a framework to track cells and groups of cells (rois) in space an and time. It allows to calculate statisitcs about lineage, shape, deformation and topology. It comes along with various tools to visulize these data. 


![alt tag](http://url/to/img.png)


How to get started?
================

We provide a preconfigured [Docker](http://docker.com/) image to run TissueMiner without any setup. Just do

     ## download some example data
     curl https://files.mpi-cbg.de/index.php/s/oEhtzFujKHUa35Z/download  | tar -zxvf -
     
     ## download the image and start the analysis
     docker pull brandl/tissue_miner
     docker run -t -i -v $(pwd)/example_movies:/movies brandl/tissue_miner /bin/bash --login -c "source /.bash_profile; cd demo_ForTissueMiner; sm all"
     
The only adjustments are to set the source directory containing your movie directories and to specify a movie of interest.


How to run locally?
================

Prerequistes. Make sure that the packages listed in [install_tm.sh](misc/install_dependencies.sh) are installed on your system.

To install TissueMiner on your machine just checkout a copy and define a TM_HOME shell variable:

    export TM_HOME="tissue_miner"

    ## download this repository
    git clone https://github.com/mpicbg-scicomp/tissue_miner.git ${TM_HOME}
        
    ## or update existing local copy with
    cd  ${TM_HOME}; git pull origin
    
    ## Install all required R packages
    ${TM_HOME}/tissue_miner/Setup.R | tee tm_setup.log
    
    ## compile the parser needed to convert TissueAnalyzer outputs into csv
    cd ${TM_HOME}/parser && make

    ## adjust your path to include all tools
    export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
    export PATH=${TM_HOME}/parser:$PATH

To execute the workflow we recommend [snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home). We provide [snakemake workflow](workflow/tm.snkmk) to ease running TissueMiner on a cluster or locally on a single computer. These are the steps to prepare your system to run TissueMiner.

    sm() {
        snakemake --snakefile ${TM_HOME}/workflow/tm.snkmk --keep-going "$@"
    }
    export -f sm
    
    cd <movie_directory>
    
    ## list all tasks
    sm -n
    
    ## process all tasks or just sepecific ones
    sm all | tee log.txt
    sm makedb | tee log.txt
    
For snakemake details see the [reference](https://bitbucket.org/johanneskoester/snakemake/wiki/Home).

Although we don not recommend it, you can run each of the tools also separately. See [simple_workflow.sh](workflow/simple_workflow.sh) for an example pipeline.

Support
=========

Please use the github ticket system to report issues or suggestions. We also welcome pull requests.