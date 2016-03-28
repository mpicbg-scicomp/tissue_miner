
# Cell contributions to tissue area change analysis

Cell contributions to tissue area changes are easily calculated from the database. It is therefore sufficient to build the database only.


### 1. Plot the rate of tissue area changes and its cellular contributions

* Copy-paste the following commands in the terminal:

```
sm make_db
cell_contributions_to_tissue_area_change_rate.R . output_analysis
```

![](cell_contributions_to_tissue_area_changes_files/figure-html/cell_contributions_to_tissue_area_changes-1.png)

[Select another analysis](tm_qs_example_data.md)


### 2. For further details
* filter by regions of interest, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#plot-the-color-coded-cell-area-pattern-in-the-whole_tissue-roi)
* compare multiple movies and ROI's, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
