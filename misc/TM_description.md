**A/** First, it requires tracked-cell data as an input. To generate such an input, we use [TissueAnalyzer](../README.md#tissueanalyzer) by [Benoit Aigouy](http://dx.doi.org/10.1016/j.cell.2010.07.042) to segment and track cells over time. **TissueMiner** contains a GUI to group cells (ROI's) in space, and a lineage browsing algorithm to follow these ROI's backward and forward in time. 

**B/** Second, **TissueMiner** consists of an **automated workflow** that stores information about:

* cell geometry
* cell topology (cell neighbor relationships)
* cell ancestry 

into a **SQLite relational database**, which is automatically queried to **quantify** and **visualize** cell dynamics during epithelium morphogenesis:

* cell state properties (position, area, anisotropy, cell packing geometry, fluorescence intensity)
* rates of cellular events (divisions, cell neighbor changes, extrusions, shape changes)
* orientation of cellular events (unit nematic description)
* rates of deformations of each type of cellular event (tensorial description)
* rates of tissue deformation (area expansion and convergence-extension) contributed by each type of cellular event (tensorial description)
* multiscale quantification and visualization using both dynamic ROI's (Lagrangian description) and fixed grids (Eulerian description): from individual cells to averages over the entire tissue

**C/** Third, **TissueMiner** also consists of a [toolkit](../README.md#documentation), which can be used in the user-friendly [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) environment, to perform *à-la-carte* analyses and visualizations of these processed data:

* time and tissue-orientation registrations for comparing multiple movies
* synchronization tools to generate a video of combined movies at constant frame rate (by under-sampling time)
* comparison of tissue and cell deformations between ROI's and between movies.
* statistics (time evolution of the distributions of cell area, anisotropy, packing, ..., bond length, vertex multiplicity, ...)
* visualization of quantified data in graphs or directly on the movie images
