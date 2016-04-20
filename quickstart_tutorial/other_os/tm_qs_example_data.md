# First Use of TissueMiner with example data

### 1. Open a Docker Quickstart Terminal

* You'll find it either in your application menu or application folder (depending on your system):

![alt tag](../../readme_screenshots/docker_toolbox_osx_quickstart_icon_nolabel.png)

### 2. Download the example data set (~100 Mb)

Just copy and paste the lines below **into this Quickstart terminal** and press Enter:
```
     ## In your terminal, type in the command below to download and extract the Demo data
     curl https://cloud.mpi-cbg.de/index.php/s/vC0VqD2Wy4A6uqu/download  | tar -zxvf -
     
     ## Go to the movie directory
     cd example_data/demo
     
     ## Store the docker command into a shell variable for the sake of simplicity in tutorials
     alias tm='docker run --rm -ti -v $(dirname $PWD):/movies -w /movies/$(basename $PWD) etournay/tissue_miner'
     
```

**A tip !** Just copy this line in your .bashrc or .bash_profile to make this command permanent:
```
alias tm='docker run --rm -ti -v $(dirname $PWD):/movies -w /movies/$(basename $PWD) etournay/tissue_miner'
```

### 3. Select the analysis you are interested in

Here, we propose some streamlined quickstart tutorials.

* [Cell area](tutorials/cell_area.md#cell-area-analysis)
* [Cell elongation](tutorials/cell_elongation.md#cell-elongation-analysis)
* [Cell packing](tutorials/cell_packing.md#cell-packing-analysis)
* [Cell lineage and divisions](tutorials/cell_lineage_and_divisions.md#cell-lineage-and-division-analysis)
* [Cell rearrangements](tutorials/cell_rearrangements.md#cell-rearrangement-analysis)
* [Cell contributions to tissue shear](tutorials/cell_contributions_to_tissue_shear.md#cell-contributions-to-tissue-shear-analysis)
* [Cell contributions to tissue area changes](tutorials/cell_contributions_to_tissue_area_changes.md#cell-contributions-to-tissue-area-change-analysis)

Here, we extend the cell area example to compare between different regions of interest.

* [Compare cell area in different ROI's](tutorials/cell_area_ROI.md#cell-area-analysis-in-multiple-rois)

Here, run an entire analysis of a single movie

```
tm sm shear_calculate topo_countt1 polygon_class tri_categorize 
tm analyze_movie.R . output_analysis
```

Here, use your own [configuration file](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/config/flywing_tm_config.R) to optimize the output rendering. Your configuration file *my_config.R* must be located in the movie repository folder. You'll find more explanation about this file on the [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#tissueminer-api-configuration)
```.
tm "export TM_CONFIG=/movies/my_config.R; sm shear_calculate topo_countt1 polygon_class tri_categorize; analyze_movie.R . output_analysis"
```


### 4. Look at the results 

* Where to find the results ? 

In the current movie folder, a **new output_analysis folder** has been created, where you'll find plots and videos created by running the tutorials. 

* Where to find the current movie folder?

In this tutorial, the example data have been downloaded into your home directory, which may be represented with a "house" icon when using a file browser. In your file browser, go to your home directory, then click on *example_data*, click on *demo*, and finally on *output_analysis*. Here, are your results !

