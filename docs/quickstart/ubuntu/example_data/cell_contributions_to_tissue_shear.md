
# Cell contributions to tissue shear analysis

Cell contributions to tissue shear (change in aspect-ratio or pure shear) are calculated from the database using a dedicated routine included in the automated workflow. Therefore, on top of building the database, we also run this routine with the command `sm shear_calculate`.


### 1. Plot the rate of tissue shear and its cellular contributions

* Copy-paste the following commands in the terminal:

```
sm shear_calculate
cell_contributions_to_tissue_shear_rate.R . output_analysis
```

![](cell_contributions_to_tissue_shear_files/figure-html/cell_contributions_to_tissue_shear_rate-1.png)

[Select another analysis](tm_qs_example_data.md)


### 2. Make videos of shear patterns

Shear pattern videos are better generated in the automated workflow. Just typing in the `sm shear_movies` command is sufficient.

* Copy-paste the following commands in the terminal:

```
sm shear_movies
```


### 3. For further details
* filter by regions of interest, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#plot-the-color-coded-cell-area-pattern-in-the-whole_tissue-roi)
* compare multiple movies and ROI's, see [Master Guide](https://mpicbg-scicomp.github.io/tissue_miner/tm_tutorial/R-tutorial.html#comparing-averaged-quantities-between-movies-and-rois)
