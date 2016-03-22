How to process a movie with TissueMiner
=====================================


First segment and track cells using TissueAnalyzer
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


Options: delimit regions of interest (ROI's) and orient the tissue along the x or y axis
=====

Fiji macros can be found in your local TissueMiner installation or [here](https://github.com/mpicbg-scicomp/tissue_miner/tree/master/fiji_macros). Just drad-and-drop the file onto Fiji and click 'run':
* draw_n_get_ROIcoord.ijm
* orient_tissue.ijm

Both macros automatically save a text file (LastFrameRoi.txt and transformation.txt, respectively). If these files are present, TissueMiner will take them into account for further precessing steps.

NB: by default, Tissue Miner will always create two ROI's: 
* "raw": corresponds to all segmented and tracked cells
* "whole_tissue": corresponds to cell lineages that remain in the field of view


Data structure
================

We invite the user to store all movies in a movie repository folder \<movie_repository\> to facilitate automated movie comparison.

Here is the **required structure** of a movie:


\<movie_repository\>/\<movie_directory\>/**Segmentation**/\<movie_directory_name\>_%03d.png

%03d represents a number padded with 3 digits. The number of digits can be modified [FAQ](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/faq.md).

All movie images are contained in the **Segmentation** folder, namely *the location images used by Tissue Analyzer* to segment and track cells.
Upon segmentation and tracking, TissueAnalyzer will generate one folder per image, in this Segmentation folder.

TissueMiner will generate additional folders and files in the \<movie_directory\> folder.


System Preparation
=======================

First, make sure to have TM_HOME shell variable pointing to the root of your TissueMiner installation as TissueMiner requires it to resolve script paths internally (see code chunck below).

Then, tweak the $PATH to include the tools (see code chunck below).

Finally, we highly recommend to add the next code block lines in your *.bashrc* or *.bash_profile* so that you don't have to export the environment variables again.

```
## adjust your path to include the tissue_miner folder
export TM_HOME="~/tissue_miner"
export

## adjust your path to include all tools for execution
export PATH=$TM_HOME/db:$TM_HOME/shear:$TM_HOME/roi:$TM_HOME/misc:$TM_HOME/movies:$TM_HOME/shear_contributions:$TM_HOME/topology:$TM_HOME/triangles:$TM_HOME/lineage:$PATH
export PATH=${TM_HOME}/parser:$PATH

## defines an alias for the snakemake command line
sm() {
    snakemake --snakefile ${TM_HOME}/workflow/tm.snkmk --keep-going "$@"
}
export -f sm

## Then apply changes
source .bashrc
```



To run TissueMiner over **your own data**
=======================

To run TissueMiner entirely over **your own data**, you'll need to go to the source directory containing your movie directories (use the *cd* command + drag and drop the folder corresponding to your *movie repository* folder onto the shell + press enter) and to replace *my_favorite_movie* with your movie directory of interest.

```
## Example for "my_favorite_movie"
cd path_to_your_movie_repository
docker run --rm -ti -v $(pwd):/movies -w /movies/my_favorite_movie brandl/tissue_miner sm all
```

For a more streamlined analysis, please visit the tutorials [here]()



Once your tracked-cell data have been processed by the TissueMiner automated workflow, you can perform a custom multiscale analysis of epithelial morphogenesis using the detailed [R-tutorial](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html) or the [Python-tutorial](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/docs/TM_tutorial_in_Python/TissueMiner_pythonTutorial-3WT_Demo.md).




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
