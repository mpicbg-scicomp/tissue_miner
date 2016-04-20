
# Cell area analysis

Cell area is stored in the database.


### 1. Make a movie of cell area pattern plotted on the tissue

* Copy-paste the following commands in the terminal:

```
tm sm make_db; tm cell_area_pattern.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_pattern-1.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 2. Plot cell area distrubution and averages
* Copy-paste the following commands in the terminal:

```
tm sm make_db; tm cell_area_graphs.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_graphs-1.png)![](cell_area_files/figure-html/cell_area_graphs-2.png)![](cell_area_files/figure-html/cell_area_graphs-3.png)![](cell_area_files/figure-html/cell_area_graphs-4.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 3. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#comparing-averaged-quantities-between-movies-and-rois)
