# First Use of TissueMiner with example data

### 1. Open a terminal

* If you don't know how to open a terminal, click [here](https://help.ubuntu.com/community/UsingTheTerminal)

### 2. Download the example data set (~100 Mb)

Just copy and paste the lines below **into your terminal** and press Enter:
```
     ## In your terminal, type in the command below to download and extract the Demo data
     curl https://cloud.mpi-cbg.de/index.php/s/EspCWSQn3k6NKLA/download  | tar -zxvf -
     
     ## Go to the movie directory
     cd example_movie/demo_ForTissueMiner
```

### 3. Select the analysis you are interested in

Here, we propose some streamlined quickstart tutorials. We systematically use the `sm` command to generate the database (if not yet present), and one dedicated R script, which takes two arguments (input movie folder and output analysis folder), to make videos or to save plots.

* [Cell area](cell_area.md)
* [Cell elongation](cell_elongation.md)
* [Cell packing](cell_packing.md)
* [Cell lineage and divisions](cell_lineage_and_divisions.md)
* [Cell rearrangements](cell_rearrangements.md)
* [Cell contributions to tissue shear](cell_contributions_to_tissue_shear.md)
* [Cell contributions to tissue area changes](cell_contributions_to_tissue_area_changes.md)


