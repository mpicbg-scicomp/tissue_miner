
# Cell contributions to tissue shear analysis

Cell contributions to tissue shear (change in aspect-ratio or pure shear) are calculated from the database using a dedicated routine included in the automated workflow. Therefore, on top of building the database, we also run this routine with the command `sm shear_calculate`.


### 1. Plot the rate of tissue shear and its cellular contributions

* Copy-paste the following commands in the terminal:

```
tm sm shear_calculate
tm cell_contributions_to_tissue_shear_rate.R . output_analysis
```

![](cell_contributions_to_tissue_shear_files/figure-html/cell_contributions_to_tissue_shear_rate-1.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)


### 2. Make videos of shear patterns

Shear pattern videos are better generated in the automated workflow. Just typing in the `sm shear_movies` command is sufficient.

* Copy-paste the following commands in the terminal:

```
tm sm shear_movies
```

* Find the videos

Open your file browser and go to the movie directory, in which you should click on the "nematics_movies" folder that contains videos of cell shear contribution pattern.

[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)


### 3. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#comparing-averaged-quantities-between-movies-and-rois)
