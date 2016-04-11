
# Cell area analysis

Cell area is stored in the database.


### 1. Make a movie of cell area pattern plotted on the tissue

* Copy-paste the following commands in the terminal:

```
$DOCKER sm make_db 
$DOCKER cell_area_pattern.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_pattern-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in) **|** 
[Try with your own data](../tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)

### 2. Plot cell area distrubution and averages
* Copy-paste the following commands in the terminal:

```
$DOCKER sm make_db 
$DOCKER cell_area_graphs.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_graphs-1.png)![](cell_area_files/figure-html/cell_area_graphs-2.png)![](cell_area_files/figure-html/cell_area_graphs-3.png)![](cell_area_files/figure-html/cell_area_graphs-4.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in) **|** 
[Try with your own data](../tm_qs_user_data.md#first-use-of-tissueminer-with-your-own-data)

### 3. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
