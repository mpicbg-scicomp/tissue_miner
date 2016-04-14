
# Cell lineage and division analysis

Cell lineages and cell divisions are stored in the database. It is therefore sufficient to build the database only.


### 1. Make a video of the cumulative pattern of cell divisions plotted on the tissue

* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_division_generation_pattern.R . output_analysis
```

![](cell_lineage_and_divisions_files/figure-html/cumulative_cell_division_pattern-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 2. Make a video of coarse-grained cell division nematics plotted on the tissue

* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_division_orientation_pattern.R . output_analysis
```

![](cell_lineage_and_divisions_files/figure-html/cell_division_nematic_pattern-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)


### 3. Plot cell division rate
* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_division_rate.R . output_analysis
```

![](cell_lineage_and_divisions_files/figure-html/cell_division_rate-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 4. Plot average cell division orientation as a function of time
* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_division_orientation.R . output_analysis
```

![](cell_lineage_and_divisions_files/figure-html/cell_division_orientation-1.png)

[How to look at the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 5. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
