
# Cell packing analysis

Cell neighbor count is easily calculated from the database. It is therefore sufficient to build the database only.


### 1. Make a video of color-coded cell neighbor number plotted on the tissue

* Copy-paste the following commands in the terminal:

```
tm sm make_db; tm cell_neighbor_number_pattern.R . output_analysis
```

![](cell_packing_files/figure-html/cell_neighbor_number_pattern-1.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 2. Plot cell neighbor count and averages
* Copy-paste the following commands in the terminal:

```
tm sm topo_countt1; tm cell_neighbor_number_graphs.R . output_analysis
```

![](cell_packing_files/figure-html/cell_neighbor_number_graphs-1.png)![](cell_packing_files/figure-html/cell_neighbor_number_graphs-2.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 3. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#comparing-averaged-quantities-between-movies-and-rois)
