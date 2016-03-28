
# Cell packing analysis

Cell neighbor count is easily calculated from the database. It is therefore sufficient to build the database only.


### 1. Make a video of color-coded cell neighbor number plotted on the tissue

* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_neighbor_number_pattern.R . output_analysis
```

![](cell_packing_files/figure-html/cell_area_pattern-1.png)

[Select another analysis](tm_qs_example_data.md)

### 2. Plot cell area distrubution and averages
* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_neighbor_number_graphs.R . output_analysis
```

![](cell_packing_files/figure-html/cell_area_graphs-1.png)![](cell_packing_files/figure-html/cell_area_graphs-2.png)

[Select another analysis](tm_qs_example_data.md)

### 3. For further details
* filter by regions of interest, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#plot-the-color-coded-cell-area-pattern-in-the-whole_tissue-roi)
* compare multiple movies and ROI's, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
