
# Cell area analysis

### 1. Make a movie of cell area pattern plotted on the tissue

We want to colour code the movie by cell area.  The `tm sm make_db` command will first build a database from the segmented data (or return "Nothing to be done" if the database is already present). The next command on the right will run the cell area analysis and make the video we want. It takes the current movie directory `.` as an input and it outputs the results in the `output_analysis` folder within the same movie directory.

* Copy-paste the following commands in the terminal:

```
tm sm make_db; tm cell_area_pattern.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_pattern-1.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 2. Plot cell area distrubution and averages

Now, instead of colour coding the cells for visualisation purposes, we want to make some graphs that demonstrate how the distributions of cell area behave in the movie. For example, we will make a graph to show how the cell area on average changes over time. The `tm sm make_db` command builds the database if not yet present (or returns "Nothing to be done" if the database is present). The next command does the plots. It takes the current movie directory `.` as an input and it outputs the results in the `output_analysis` folder within the same movie directory.

* Copy-paste the following commands in the terminal:

```
tm sm make_db; tm cell_area_graphs.R . output_analysis
```

![](cell_area_files/figure-html/cell_area_graphs-1.png)![](cell_area_files/figure-html/cell_area_graphs-2.png)![](cell_area_files/figure-html/cell_area_graphs-3.png)![](cell_area_files/figure-html/cell_area_graphs-4.png)

[Where to find the results ?](../tm_qs_example_data.md#4-look-at-the-results) **|** 
[Back to tutorial list](../tm_qs_example_data.md#3-select-the-analysis-you-are-interested-in)

### 3. For further details

* compare multiple movies and ROI's, see [TM R User Manual](https://mpicbg-scicomp.github.io/tissue_miner/user_manual/TM_R-UserManual.html#comparing-averaged-quantities-between-movies-and-rois)
