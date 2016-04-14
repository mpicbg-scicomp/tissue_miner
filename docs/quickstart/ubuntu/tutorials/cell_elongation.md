
# Cell elongation analysis

Cell elongation is stored in the database. It is therefore sufficient to build the database only.


### 1. Make a movie of cell elongation magnitude pattern plotted on the tissue

* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_elongation_magnitude_pattern.R . output_analysis
```

![](cell_elongation_files/figure-html/cell_elongation_magnitude_pattern-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)


### 2. Make a movie of coarse-grained cell elongation nematics plotted on the tissue
* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_elongation_nematic_cg_pattern.R . output_analysis
```

![](cell_elongation_files/figure-html/cg_cell_elongation_nematic_pattern-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 3. Plot cell elongation distrubution and averages
* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_elongation_graphs.R . output_analysis
```

![](cell_elongation_files/figure-html/cell_elongation_graphs-1.png)![](cell_elongation_files/figure-html/cell_elongation_graphs-2.png)![](cell_elongation_files/figure-html/cell_elongation_graphs-3.png)![](cell_elongation_files/figure-html/cell_elongation_graphs-4.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 4. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
