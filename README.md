
About
=================

**TissueMiner** is a framework to track cells and groups of cells (rois) in space an and time. It allows to calculate statistics about lineage, shape, deformation and topology. It comes along with various tools to visualize these data. 


![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/stripes_0.jpg)
![alt tag](https://raw.githubusercontent.com/mpicbg-scicomp/tissue_miner/master/docs/readme_screenshots/veins_0.jpg)


How to get started?
================

We provide a preconfigured [docker image](https://registry.hub.docker.com/u/brandl/tissue_miner/) to run TissueMiner without any setup. Just do

     ## download some example data
     curl https://files.mpi-cbg.de/index.php/s/oEhtzFujKHUa35Z/download  | tar -zxvf -
     
     ## download the image and start the analysis
     docker pull brandl/tissue_miner
     docker run -t -i -v $(pwd)/example_movies:/movies brandl/tissue_miner /bin/bash --login -c "source /.bash_profile; cd demo_ForTissueMiner; sm all"
     
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
    cd  ${TM_HOME}; git pull origin
    
    ## Install all required R packages
    ${TM_HOME}/tissue_miner/Setup.R | tee tm_setup.log
    
    ## compile the parser needed to convert TissueAnalyzer outputs into csv
    cd ${TM_HOME}/parser && make

    ## adjust your path to include all tools
    export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
    export PATH=${TM_HOME}/parser:$PATH
    
Don't forget to define a TM_HOME shell variable, pointing to the root of your TissueMiner installation, since it will require it to resolve script paths internally.

Then you can run TissueMiner to analyze your movie provided you movie folder contains a "Segmentation" folder in which to find the original images as well as the TissueAnalyzer outputs.

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
    sm makedb | tee log.txt
    
    ## Export statistics and execution state graph visualization
    sm --dag | dot -Tpdf > dag_tbd.pdf
    sm -D > sm_execution_state.txt
    
For snakemake details see the its [reference](https://bitbucket.org/johanneskoester/snakemake/wiki/Home).

Although we do not recommend it, you can also run each of the TissueMiner tools separately. See [simple_workflow.sh](workflow/simple_workflow.sh) for an example pipeline.

Support
=========

Please use the github ticket system to report issues or suggestions. We also welcome pull requests.