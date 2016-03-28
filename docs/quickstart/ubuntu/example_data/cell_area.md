
# Cell area analysis

Cell area is stored in the database.


### 1. Make a movie of cell area pattern plotted on the tissue

* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_area_pattern.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_pattern-1.png)

[Select another analysis](tm_qs_example_data.md)

### 2. Plot cell area distrubution and averages
* Copy-paste the following commands in the terminal:

```
sm make_db 
cell_area_graphs.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_graphs-1.png)![](cell_area_files/figure-html/cell_area_graphs-2.png)![](cell_area_files/figure-html/cell_area_graphs-3.png)![](cell_area_files/figure-html/cell_area_graphs-4.png)

[Select another analysis](tm_qs_example_data.md)

### 3. For further details

* compare multiple movies and ROI's, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
