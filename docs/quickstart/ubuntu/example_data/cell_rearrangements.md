
# Cell rearrangement analysis

Cell neighbor changes are calculated from the database using a dedicated routine included in the automated workflow. Therefore, on top of building the database, we also run this routine with the command `sm topo_countt1`


### 1. Make a video of the cell rearrangement coarse-grained pattern plotted on the tissue

* Copy-paste the following commands in the terminal:

```
sm topo_countt1 
cell_neighbor_change_orientation_pattern.R . output_analysis
```

![](cell_rearrangements_files/figure-html/cell_neighbor_change_orientation_pattern-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in) **|** 
[Try with your own data](../tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)


### 2. Plot cell neighbor change rate
* Copy-paste the following commands in the terminal:

```
sm topo_countt1 
cell_neighbor_change_rate.R . output_analysis
```

![](cell_rearrangements_files/figure-html/cell_neighbor_change_rate-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in) **|** 
[Try with your own data](../tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)

### 3. Plot average cell neighbor change orientation as a function of time
* Copy-paste the following commands in the terminal:

```
sm topo_countt1 
cell_neighbor_change_orientation.R . output_analysis
```

![](cell_rearrangements_files/figure-html/cell_neighbor_change_orientation-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in) **|** 
[Try with your own data](../tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)

### 4. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
