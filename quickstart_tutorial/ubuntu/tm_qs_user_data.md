# First Use of TissueMiner with your own data

### 1. First, organize your movie images as follows

* Create a **movie_repository folder** where you store the data from all your movies
    + Inside this folder, create one **movie folder** per movie
        + inside each movie folder, create a **Segmentation folder** to store movie images
            + the movie should be stored as a series of tif or png files, one per timepoint, with the name of each file composed of a sample name + underscore + frame number 

The data organization is summarized as follow:

`<movie_repository>/<movie_directory>/Segmentation/<movie_directory_name>_%03d.png`


where %03d represents a frame number padded with 3 digits. The number of digits can be modified, see [FAQ](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/faq.md).


Here, is an example: movieSegmentation/demo/Segmentation/demo_000.png

![alt tag](https://github.com/mpicbg-scicomp/tissue_miner/blob/gh-pages/readme_screenshots/data_organization.png)

* The image type should be **8-bit** or **RGB**.
* The image format should be **png** or **tiff**


### 2. Segment and track cells using TissueAnalyzer

**Tissue Analyzer** (copyright Aigouy 2016) ships with its own licence
(see license_TA.txt bundled in the software). Tissue Analyzer should not be modified or
reverse engineered.  Tissue Analyzer should always be distributed
bundled with **TissueMiner** and not alone.

You can install the latest version of Tissue Analyzer (formerly known
as Packing Analyzer) within FIJI (http://fiji.sc/Fiji). To do so:

* get a fresh FIJI installation (including JDK8)
* launch FIJI
* open the "Help" menu
* click on "Update..."
* click on "Manage update sites"
* click on "Add" and enter "http://sites.imagej.net/TA/" in the URL field (and anything you like in the "Name" field)
* click "Close" and FIJI should offer you to install TA.

Once the installation is complete, restart FIJI, open the "Plugins" menu and click on "Tissue Analyzer"

For a quick start guide, click [here](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/docs/TAdoc.pdf)


### 3. Optionally, define regions of interest (ROI's) and orient the tissue along the x or y axis
To this purpose, we provide two [FIJI](http://fiji.sc/) programs ***draw_n_get_ROIcoord.ijm*** and ***orient_tissue.ijm***.

#### 3.1 Define ROI's:
* launch FIJI
* go to the ***fiji_macros*** folder located in your TissueMiner installation folder
* drag-and-drop the ***draw_n_get_ROIcoord.ijm*** file into FIJI
* a script editor opens automatically
* click **RUN** and define your ROI's   

By default, TissueMiner always creates two ROI's: 

* "raw": corresponds to all segmented and tracked cells
* "whole_tissue": corresponds to cell lineages that remain in the field of view



#### 3.2 Orient the tissue
* launch FIJI
* go to the ***fiji_macros*** folder located in your TissueMiner installation folder
* drap-and-drop the ***orient_tissue.ijm*** file into FIJI
* a script editor opens automatically
* click **RUN** and define the tissue axis to be aligned on x or y. 

Both programs automatically save a text file (UserFrameRoi.txt and transformation.txt, respectively). If these files are present, TissueMiner will take them into account for further precessing steps.


### 4. Open a terminal

* If you don't know how to open a terminal, click [here](https://help.ubuntu.com/community/UsingTheTerminal)

### 5. Go to the movie directory

* Browse your data to find the **movie folder**
* **In your terminal** type in ***cd*** and press **space**
* and drag-and-drop the **movie folder** into the terminal

Example:
`cd /home/tissueminer/movie_repository/movie_1` 

### 6. Select the analysis you are interested in
Here, we propose some streamlined quickstart tutorials.

* [Cell area](tutorials/cell_area.md#cell-area-analysis)
* [Cell elongation](tutorials/cell_elongation.md#cell-elongation-analysis)
* [Cell packing](tutorials/cell_packing.md#cell-packing-analysis)
* [Cell lineage and divisions](tutorials/cell_lineage_and_divisions.md#cell-lineage-and-division-analysis)
* [Cell rearrangements](tutorials/cell_rearrangements.md#cell-rearrangement-analysis)
* [Cell contributions to tissue shear](tutorials/cell_contributions_to_tissue_shear.md#cell-contributions-to-tissue-shear-analysis)
* [Cell contributions to tissue area changes](tutorials/cell_contributions_to_tissue_area_changes.md#cell-contributions-to-tissue-area-change-analysis)

All tutorials can be used to compare between **different regions of interest**. Here, we give an example with cell area

* [Compare cell area in different ROI's](tutorials/cell_area_ROI.md#cell-area-analysis-in-multiple-rois)

Here, run an entire analysis of a single movie

```
sm shear_calculate topo_countt1 polygon_class tri_categorize 
analyze_movie.R . output_analysis
```

Here, use your own [configuration file](https://github.com/mpicbg-scicomp/tissue_miner/blob/master/config/flywing_tm_config.R) to optimize the output rendering. Your configuration file *my_config.R* must be located in the movie repository folder. You'll find more explanation about this file in the [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#tissueminer-api-configuration).
```
export TM_CONFIG=$(dirname $PWD)/my_config.R; sm shear_calculate topo_countt1 polygon_class tri_categorize; analyze_movie.R . output_analysis
```


### 7. Look at the results 

* Where to find the results ? 

In the current movie folder, a **new output_analysis folder** has been created, where you'll find plots and videos created by running the tutorials. 

* Where to find the current movie folder?

In this tutorial, the example data have been downloaded into your home directory, which may be represented with a "house" icon when using a file browser. In your file browser, go to your home directory, then click on *example_data*, click on *demo*, and finally on *output_analysis*. Here, are your results !
