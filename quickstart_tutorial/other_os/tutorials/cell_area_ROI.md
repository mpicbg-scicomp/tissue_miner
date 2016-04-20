
# Cell area analysis in multiple ROI's

Cell area is stored in the database, but ROI's are not. ROI's are identified by their name. Here are the available ROI's for the example data: 

* raw 
* whole_tissue
* cell_patch

### 1. Make a movie of cell area pattern plotted on the tissue for the *cell_patch* ROI

* Copy-paste the following commands in the terminal:

```
tm "sm roi_tracking; cell_area_pattern.R . output_analysis 'cell_patch'"
```

![](cell_area_ROI_files/figure-html/cell_area_pattern-1.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 2. Plot cell area distrubution and averages in each ROI
* Copy-paste the following commands in the terminal:

```
tm "sm make_db; cell_area_graphs.R . output_analysis 'raw whole_tissue cell_patch'"
```

![](cell_area_ROI_files/figure-html/cell_area_graphs-1.png)![](cell_area_ROI_files/figure-html/cell_area_graphs-2.png)![](cell_area_ROI_files/figure-html/cell_area_graphs-3.png)![](cell_area_ROI_files/figure-html/cell_area_graphs-4.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 3. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#comparing-averaged-quantities-between-movies-and-rois)
