
# TissueMiner: Python tutorial
This tutorial is intended to demonstrate how Python can be used for visualisation and quantification of the data produced by the TissueMiner workflow. It can be run as Ipython Notebook.

Prerequisites:
    - time-lapse movies have been processed using TissueMiner workflow 
    - python libraries: numpy, matplotlib, pandas, ipython, sqlite3, rpy2
    - provided tissue_miner.py and tisue_miner_tools.py files are stored in the same folder -> no other installation is necessary

Tutorial is centered around Pandas and Matplotlib libraries. Although we provide several simple examples, basic knowledge of python and these two libraries is recommended. For a 10 minutes introduction to Pandas see: http://pandas.pydata.org/pandas-docs/stable/10min.html . Beginner guide to Matplotlib can be found here: http://matplotlib.org/users/beginner.html


Tutorial is organized in three parts:
    - Introduction
        - How to use Pandas and Matplotlib libraries for basic analysis of TissueMiner data?
    - Visualisation
        - How to visualise different cell properties on original images in cellular resolution?
    - Comparing movies and ROIs
        - How to compare data between different movies and/or regions of interest (ROIs)?

# Introduction


```python
# Set path to tissue_miner.py and tissue_miner_tools.py
libPath= 'lib/'

import sys
sys.path.insert(0, libPath)

import tissue_miner_tools as tml          
import tissue_miner as tm
```


```python
import os 

# Set path to the movie databases
movieDatabaseDir='/data/biophys/mpopovic/temp/MovieRepository_DB/'

# Set relative path to the ROI file
movieRoiPath= 'roi_bt/' 
movieRoiFile= 'lgRoiSmoothed'

# Set path to output folder (optional)
outDatabaseDir= '/home/mpopovic/Dropbox/resource_paper/tutorial_output'

if not os.path.exists(outDatabaseDir):
    os.makedirs(outDatabaseDir)
```


```python
# matplotlib magic function to enable plotting in the Notebook
%matplotlib inline  

# Import (some of) prerequisite libraries, others are used by tissue_miner and tissue_miner_tools
import os as os

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import collections as mc

from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pandas as pd
import pandas.io.sql as psql
```

## Database and supplementary data
### Movie class
Movie class is the basic object for manipulating the data produced by TissueMiner. It is implemented in tissue_miner.py and can be simply extended by user. 
Each movie will have a dedicated instance of Movie class which is supplied with methods for loading and performing basic operations on the corresponding movie data. 

Provided version of the Movie class is provided with methods for loading SQL database tables as Pandas DataFrames. On initialization it loads cells and frames tables as well as region of interests (ROI) table if provided and reads time shifts provided in file 'PupalWingMovies.csv' which are used to synchronize movies in time.

Note that movie class will by default create pickle (.pkl) format version of both SQL tables and extra tables (see below) on the first run. This is done to reduce loading time of the data in future runs. 


```python
name= 'WT_2'

# Initialization
movie= tm.Movie(name, path= movieDatabaseDir, ROI_path= movieRoiPath, ROI_filename= movieRoiFile) 
```

    Loading table frames from database WT_2...
    Loading table cells from pickle file...
    Loading roiBT ...


### Loading the database tables
Movie method load_DB_table(table_name) loads table named table_name from the database. Table is loaded as Pandas DataFrame object and stored in dictionary movie.DB_table.


```python
movie.load_DB_table('cells')

movie.DB_table['cells'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>cell_id</th>
      <th>center_x</th>
      <th>center_y</th>
      <th>area</th>
      <th>elong_xx</th>
      <th>elong_xy</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>10001</td>
      <td>195.161416</td>
      <td>1120.430486</td>
      <td>513.0</td>
      <td>0.150421</td>
      <td>-0.002361</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>10001</td>
      <td>187.481634</td>
      <td>1118.081728</td>
      <td>506.0</td>
      <td>0.125004</td>
      <td>-0.003316</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>10001</td>
      <td>178.082442</td>
      <td>1113.230519</td>
      <td>492.5</td>
      <td>0.155494</td>
      <td>0.040019</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>10001</td>
      <td>172.813431</td>
      <td>1109.035565</td>
      <td>538.5</td>
      <td>0.142844</td>
      <td>0.094702</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>10001</td>
      <td>167.999972</td>
      <td>1106.556618</td>
      <td>488.0</td>
      <td>0.198969</td>
      <td>0.125119</td>
    </tr>
  </tbody>
</table>
</div>



### Loading the supplement data
In addition to the SQL database, extra tables are provided for each time-lapse movie. Each set of supplment data has a dedicated method for loading. In this tutorial we will use:
    - cellshapes.RData - cell contours by using anticlockwisely ordered cell vertices
    - ./roi_bt/lgRoiSmoothed.RData - user-defined regions of interest
    - ./shear_contrib/<ROI_name>/avgDeformTensorsWide.tsv - precalculated deformation rates of triangles and tissue for each region of interest

Note that some of these extra tables are provided in .RData format. On the first run these are converted to pickle format which might take some time.


```python
# For example, loading cellshapes table:
movie.load_cellshapes()

movie.cellshapes.head()
```

    Loading cellshapes of WT_2...





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>cell_id</th>
      <th>x_pos</th>
      <th>y_pos</th>
      <th>bond_order</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>10001</td>
      <td>179.102958</td>
      <td>1118.318960</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0</td>
      <td>10001</td>
      <td>178.519784</td>
      <td>1125.365939</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0</td>
      <td>10001</td>
      <td>200.016491</td>
      <td>1133.104904</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0</td>
      <td>10001</td>
      <td>210.582070</td>
      <td>1125.464974</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0</td>
      <td>10001</td>
      <td>205.758975</td>
      <td>1111.726788</td>
      <td>5</td>
    </tr>
  </tbody>
</table>
</div>



## Manipulation of large data sets using Pandas library
Here we provide several simple examples of data manipulation using Pandas 


```python
# We first convert cell area of each cell from pixel squared to micrometer squared by creating a new column 
# named 'realArea' to DB_table['cells'] DataFrame by multiplying 'area' column with 
# conversion factor 1pixel = 0.207 micrometers
movie.DB_table['cells']['realArea']= movie.DB_table['cells']['area']*0.207**2

movie.DB_table['cells'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>cell_id</th>
      <th>center_x</th>
      <th>center_y</th>
      <th>area</th>
      <th>elong_xx</th>
      <th>elong_xy</th>
      <th>realArea</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>10001</td>
      <td>195.161416</td>
      <td>1120.430486</td>
      <td>513.0</td>
      <td>0.150421</td>
      <td>-0.002361</td>
      <td>21.981537</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>10001</td>
      <td>187.481634</td>
      <td>1118.081728</td>
      <td>506.0</td>
      <td>0.125004</td>
      <td>-0.003316</td>
      <td>21.681594</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>10001</td>
      <td>178.082442</td>
      <td>1113.230519</td>
      <td>492.5</td>
      <td>0.155494</td>
      <td>0.040019</td>
      <td>21.103132</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>10001</td>
      <td>172.813431</td>
      <td>1109.035565</td>
      <td>538.5</td>
      <td>0.142844</td>
      <td>0.094702</td>
      <td>23.074186</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>10001</td>
      <td>167.999972</td>
      <td>1106.556618</td>
      <td>488.0</td>
      <td>0.198969</td>
      <td>0.125119</td>
      <td>20.910312</td>
    </tr>
  </tbody>
</table>
</div>




```python
# We then calculate average cell area in each frame. To this end we select 'frame' and 'area columns of DB_table['cells'] 
# Dataframe and we group rows by values of 'frame' column. On each of the groups mean value of area is calculated 
movie.average_cell_area= movie.DB_table['cells'][['frame', 'realArea']].groupby('frame').agg(np.mean).reset_index().sort(columns= 'frame')

movie.average_cell_area.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>realArea</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>25.231777</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>24.815621</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>24.307623</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>23.862201</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>23.441276</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Table 'frames' provides time elapsed between since the first frame of the movie in seconds. It is loaded by default.

movie.DB_table['frames'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>time_sec</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>280</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>561</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>841</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>1124</td>
    </tr>
  </tbody>
</table>
</div>




```python
# We can now add time information to the average_cell_area DataFrame by merging it with movie.DB_table['frames'] DataFrame.
# Merging is done on column 'frame' present in both DataFrames.
movie.average_cell_area= movie.average_cell_area.merge(movie.DB_table['frames'], on= 'frame')

movie.average_cell_area.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>realArea</th>
      <th>time_sec</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>25.231777</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>24.815621</td>
      <td>280</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>24.307623</td>
      <td>561</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>23.862201</td>
      <td>841</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>23.441276</td>
      <td>1124</td>
    </tr>
  </tbody>
</table>
</div>




```python
# New column 'time' is defined to store time information in hours
movie.average_cell_area['time']= movie.average_cell_area['time_sec']/3600.

movie.average_cell_area.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>realArea</th>
      <th>time_sec</th>
      <th>time</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>25.231777</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>24.815621</td>
      <td>280</td>
      <td>0.077778</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>24.307623</td>
      <td>561</td>
      <td>0.155833</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>23.862201</td>
      <td>841</td>
      <td>0.233611</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>23.441276</td>
      <td>1124</td>
      <td>0.312222</td>
    </tr>
  </tbody>
</table>
</div>




```python
# We can now use DataFrame movie.average_cell_area to visualize average cell area in time using matplotlib.pyplot 
plt.plot(movie.average_cell_area['time'], movie.average_cell_area['realArea'])

plt.ylabel(r'area $[\mu m^2]$')
plt.xlabel(r'time $[h]$')
plt.title('Average area as a function of time')
plt.ylim(13, 26)
plt.grid()
plt.savefig(outDatabaseDir+'averageAreaInTime.png')
```


![png](output_17_0.png)


# Visualization
Here we demonstrate how to load and visualise cell data on the original images.
We first show all the details on a simple example and then we use Movie.plot_frame_cells() method to visualize cell area, cell elongation, cell packnig and division patterns. All examples on are done for a single frame.


```python
# We first demonstrate how to visualize cell outlines using extra table cellshapes
frame= 100

# To visualise cells we subset cellshapes table (loaded in Introduction) by the frame
movie.frame_cellshapes= tml.subset_dataframe(movie.cellshapes, ['frame'], [frame])

# Construct a list of polygon vertices for each cell
movie.frame_polygons= movie.frame_cellshapes.groupby('cell_id').apply(lambda x: zip(x['x_pos'].values, x['y_pos'].values)).reset_index().rename(columns= {0: 'plot_vertices'})

# We use matplotlib.collections to visualise a large number of polygons at once
plt.figure(figsize= (6, 10))    # define figure size
movie.show_image(frame)          # show original image of the frame     
plt.gca().autoscale_view()      
plt.gca().set_aspect('equal')   
coll= mc.PolyCollection(movie.frame_polygons['plot_vertices'].values, lw= .5)   # create a collection of polygons
coll.set_edgecolor('lightgreen')
coll.set_facecolor('none')
plt.gca().add_collection(coll)        # add collection to the figure
plt.xlim(0, 900) 
plt.ylim(300, 1800) 
plt.savefig(outDatabaseDir+'cell_outlines_reversed.png')
```

    /usr/lib64/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
      if self._edgecolors == str('face'):



![png](output_19_1.png)



```python
# Image is revesred because coordinate system of the image has y-axis pointing down - we have to reverse y-axis!

frame= 100

movie.frame_cellshapes= tml.subset_dataframe(movie.cellshapes, ['frame'], [frame])

movie.frame_polygons= movie.frame_cellshapes.groupby('cell_id').apply(lambda x: zip(x['x_pos'].values, x['y_pos'].values)).reset_index().rename(columns= {0: 'plot_vertices'})

plt.figure(figsize= (6, 10))
movie.show_image(frame)                                                      
plt.gca().autoscale_view()                                                  
plt.gca().set_aspect('equal')                                                
coll= mc.PolyCollection(movie.frame_polygons['plot_vertices'].values, lw= .5)         
coll.set_edgecolor('lightgreen')
coll.set_facecolor('none')
plt.gca().add_collection(coll)                                              
plt.xlim(0, 900)                                                           
plt.ylim(300, 1800)                                                                                                              
plt.axis('off')      
plt.gca().invert_yaxis()  # Inverting the y-axis of the plot
plt.savefig(outDatabaseDir+'cell_outlines.png')
 
```


![png](output_20_0.png)


## Cell area


```python
# Preparind the data in frame 70 for visualization.
frame= 70

# To visualise cells we subset cellshapes table (loaded in Introduction) by the frame
movie.frame_cellshapes= tml.subset_dataframe(movie.cellshapes, ['frame'], [frame])

# Construct a list of polygon vertices for each cell
movie.frame_polygons= movie.frame_cellshapes.groupby('cell_id').apply(lambda x: zip(x['x_pos'].values, x['y_pos'].values)).reset_index().rename(columns= {0: 'plot_vertices'})

# We subset movie.DB_table['cells'] DataFrame to select cells in a given frame
movie.frame_cells= tml.subset_dataframe(movie.DB_table['cells'], ['frame'], [frame])

# Add the cell area information to the frame_polygons table
movie.frame_polygons= movie.frame_polygons.merge(movie.frame_cells[['cell_id', 'realArea']], on= 'cell_id')
```


```python
# To visualise cell area we use method Movie.plot_frame_cells()
movie.plot_frame_cells(frame, movie.frame_polygons, title= 'cell area $[\mu m^2]$', color_column= 'realArea', c_min= 0, c_max= 50, color_map= cm.gist_rainbow)
plt.savefig(outDatabaseDir+'cell_area_frame_' + tml.fill_zeros(str(frame), 3) + '.png')
```


![png](output_23_0.png)


## Cell elongation


```python
# Preparind the data in frame 70 for visualization.
frame= 70

movie.frame_cellshapes= tml.subset_dataframe(movie.cellshapes, ['frame'], [frame])

movie.frame_polygons= movie.frame_cellshapes.groupby('cell_id').apply(lambda x: zip(x['x_pos'].values, x['y_pos'].values)).reset_index().rename(columns= {0: 'plot_vertices'})

movie.frame_cells= tml.subset_dataframe(movie.DB_table['cells'], ['frame'], [frame])

movie.frame_polygons= movie.frame_polygons.merge(movie.frame_cells[['cell_id', 'elong_xx']], on= 'cell_id')
```


```python
movie.plot_frame_cells(frame, movie.frame_polygons, title= 'PD cell elongation', color_column= 'elong_xx', c_min= 0, c_max= .7, color_map= cm.gist_rainbow)
plt.savefig(outDatabaseDir+'cell_elong_frame_' + tml.fill_zeros(str(frame), 3) + '.png')
```


![png](output_26_0.png)


## Cell packing 


```python
# To determine number of neighbor for each cell we need to load 'directed bonds' table and store it 
# in movie.DB_table['directed_bonds'] DataFrame
movie.load_DB_table('directed_bonds')

movie.DB_table['directed_bonds'].head()
```

    Loading table directed_bonds from pickle file...





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>frame</th>
      <th>cell_id</th>
      <th>dbond_id</th>
      <th>conj_dbond_id</th>
      <th>bond_id</th>
      <th>vertex_id</th>
      <th>left_dbond_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>10001</td>
      <td>1130</td>
      <td>1129</td>
      <td>567</td>
      <td>393</td>
      <td>1312</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0</td>
      <td>10001</td>
      <td>1236</td>
      <td>1235</td>
      <td>620</td>
      <td>427</td>
      <td>1130</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0</td>
      <td>10001</td>
      <td>1312</td>
      <td>1311</td>
      <td>658</td>
      <td>392</td>
      <td>1439</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0</td>
      <td>10001</td>
      <td>1439</td>
      <td>1441</td>
      <td>722</td>
      <td>452</td>
      <td>1442</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0</td>
      <td>10001</td>
      <td>1442</td>
      <td>1443</td>
      <td>723</td>
      <td>493</td>
      <td>1444</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Preparing the data in frame 70 for visualization.
frame= 70

# We subset movie.DB_table['directed_bonds'] DataFrame to select bonds in a given frame
movie.frame_dbonds= tml.subset_dataframe(movie.DB_table['directed_bonds'], ['frame'], [frame])
```


```python
# Number of neighbors is equal to number of bonds of each cell. Therefore we group rows of movie.frame_dbonds
# DataFrame by cell_id and we assing to each cell_id number of element in the corresponding group
movie.frame_cell_packing= movie.frame_dbonds[['cell_id', 'dbond_id']].groupby('cell_id').agg(len).reset_index().rename(columns= {'dbond_id': 'nr_neighbors'})

movie.frame_cell_packing.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cell_id</th>
      <th>nr_neighbors</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>10002</td>
      <td>5</td>
    </tr>
    <tr>
      <th>1</th>
      <td>10003</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10011</td>
      <td>5</td>
    </tr>
    <tr>
      <th>3</th>
      <td>10022</td>
      <td>7</td>
    </tr>
    <tr>
      <th>4</th>
      <td>10033</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>




```python
movie.frame_polygons= movie.frame_cellshapes.groupby('cell_id').apply(lambda x: zip(x['x_pos'].values, x['y_pos'].values)).reset_index().rename(columns= {0: 'plot_vertices'})

# As before, we add information about number of neighbors to movie.frame_polygons DataFrame 
movie.frame_polygons= movie.frame_polygons.merge(movie.frame_cell_packing, on= 'cell_id')
```


```python
movie.plot_frame_cells(frame, movie.frame_polygons, title= 'cell packing', color_column= 'nr_neighbors', c_min= 3., c_max= 8., color_map= cm.gist_rainbow)
plt.savefig(outDatabaseDir+'cell_packnig_frame_' + tml.fill_zeros(str(frame), 3) + '.png')
```


![png](output_32_0.png)


## Cell division patterns


```python
# To visualize division patterns we will show to which generation each cell belongs. This information is  
# available in 'cell_histories' table of SQL database.

movie.load_DB_table('cell_histories')

movie.DB_table['cell_histories'].head()
```

    Loading table cell_histories from pickle file...





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>cell_id</th>
      <th>tissue_analyzer_group_id</th>
      <th>first_occ</th>
      <th>last_occ</th>
      <th>left_daughter_cell_id</th>
      <th>right_daughter_cell_id</th>
      <th>appears_by</th>
      <th>disappears_by</th>
      <th>lineage_group</th>
      <th>generation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>10001</td>
      <td>1265</td>
      <td>0</td>
      <td>20</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Unclassified</td>
      <td>Apoptosis</td>
      <td>lg_2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>10002</td>
      <td>28170</td>
      <td>0</td>
      <td>177</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Unclassified</td>
      <td>MovesOutOfMask</td>
      <td>lg_3</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10003</td>
      <td>52814</td>
      <td>0</td>
      <td>106</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Unclassified</td>
      <td>MovesOutOfMask</td>
      <td>lg_4</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>10004</td>
      <td>64282</td>
      <td>0</td>
      <td>3</td>
      <td>11124</td>
      <td>11125</td>
      <td>Unclassified</td>
      <td>Division</td>
      <td>lg_5</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>10005</td>
      <td>78460</td>
      <td>0</td>
      <td>6</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Unclassified</td>
      <td>MovesOutOfMask</td>
      <td>lg_6</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Combine generatino information with frame_cells table

# In this example we show late because it shows cumulative division pattern.
frame= 170

movie.frame_cellshapes= tml.subset_dataframe(movie.cellshapes, ['frame'], [frame])

movie.frame_polygons= movie.frame_cellshapes.groupby('cell_id').apply(lambda x: zip(x['x_pos'].values, x['y_pos'].values)).reset_index().rename(columns= {0: 'plot_vertices'})

movie.frame_polygons= movie.frame_polygons.merge(movie.DB_table['cell_histories'][['cell_id', 'generation']], on= 'cell_id')
```


```python
movie.plot_frame_cells(frame, movie.frame_polygons, title= 'cell generation', color_column= 'generation', c_min= 0., c_max= 3., n_ticks= 3,color_map= cm.gist_rainbow)
plt.savefig(outDatabaseDir+'cell_division_pattern_frame_' + tml.fill_zeros(str(frame), 3) + '.png')
```


![png](output_36_0.png)


# Comparing average quantities among movies and ROIs
Data produced by TissueMiner is easy to compare among different movies and regions of interest. Here we first demonstrate the comparison of total cell number in three provided movies and then we compare average cell area and elongation between the three movies in different regions of interest.


```python
# Define a list of movies to compare

movie_list= ['WT_1',
             'WT_2',
             'WT_3']

# Define a list of ROIs to analyse
ROI_list= ['whole_tissue', 'distL3', 
           'distInterL2-L3', 'distInterL3-L4']
```


```python
movies= {}  # Dictionary to store instances of Movie class for different movies - contains all the data we will use.

for name in movie_list:
    print name
    movies[name]= tm.Movie(name, path= movieDatabaseDir, ROI_path= movieRoiPath, ROI_filename= movieRoiFile)
    
```

    WT_1
    Loading table frames from database WT_1...
    Loading table cells from pickle file...
    Loading roiBT ...
    WT_2
    Loading table frames from database WT_2...
    Loading table cells from pickle file...
    Loading roiBT ...
    WT_3
    Loading table frames from database WT_3...
    Loading table cells from pickle file...
    Loading roiBT ...


## Cell count for different movies


```python
# We determine number of cells in 'blade' and 'whole_tissue' regions of interest using method Movie.cell_number(region).
for name in movie_list:
    movies[name].blade_nr_cells= tml.cells_number(movies[name].region_cells('distInterL3-L4'))
    movies[name].whole_tissue_nr_cells= tml.cells_number(movies[name].region_cells('whole_tissue'))

plt.figure(figsize= (12, 6))    
plt.subplot(1, 2, 1)
for name in movie_list:
    plt.plot(movies[name].time[:movies[name].NFrames], movies[name].blade_nr_cells, label= name)
plt.grid()
plt.title('distInterL3-L4')
plt.xlabel('Time [hAPF]')
plt.ylabel('Cell count')
plt.subplot(1, 2, 2)
for name in movie_list:
    plt.plot(movies[name].time[:movies[name].NFrames], movies[name].whole_tissue_nr_cells, label= name)
plt.grid()
plt.title('whole_tissue')
plt.xlabel('Time [hAPF]')
plt.ylabel('Cell count')
plt.legend(loc= 'best')
plt.tight_layout()
plt.savefig(outDatabaseDir+'cell_count_comparison.png')
```


![png](output_41_0.png)


## Average cell area in different ROIs and movies


```python
# Calculate average cell area for each ROI of each movie using region_cells_area_avg from tissue_miner_tools.py.
# Conversion factor 0.207 from pixels to micrometers is included.

for name in movie_list:
    movies[name].average_cell_area= {}

for region in ROI_list:
    for name in movie_list:
        movies[name].average_cell_area[region]= tml.cells_average_area(movies[name].region_cells(region))*.207**2
        
```


```python
plt.figure(figsize= (16, 10))
for region in ROI_list:
    plt.subplot(3, 4, 1 + ROI_list.index(region))
    for name in movie_list:
        plt.plot(movies[name].time[:movies[name].NFrames], movies[name].average_cell_area[region], label= name)
    plt.ylim(5, 40)
    plt.xlabel('Time [hAPF]')
    plt.ylabel('cell area $[\mu m^2]$')
    plt.title(region)
    plt.grid()
plt.legend(loc= 'best')
plt.tight_layout()
plt.savefig(outDatabaseDir+'cell_average_area_ROI_comparison.png')
```


![png](output_44_0.png)


## Norm of  average cell elongation in different ROIs and movies


```python
# Calculate average cell elongation for each ROI of each movie using region_cel_shape_avg() function from tissue_miner_tools.py

for name in movie_list:
    movies[name].average_cell_elong_xx= {}
    movies[name].average_cell_elong_xy= {}

for region in ROI_list:
    for name in movie_list:
        (movies[name].average_cell_elong_xx[region], 
         movies[name].average_cell_elong_xy[region])= tml.cells_average_elongation(movies[name].region_cells(region))
```


```python
plt.figure(figsize= (16, 10))
for region in ROI_list:
    plt.subplot(3, 4, 1 + ROI_list.index(region))
    for name in movie_list:
        plt.plot(movies[name].time[:movies[name].NFrames], 
                 np.sqrt(movies[name].average_cell_elong_xx[region]**2 +
                         movies[name].average_cell_elong_xy[region]**2), label= name)
        plt.ylim(0, 0.32)
    plt.title(region)
    plt.xlabel('Time [hAPF]')
    plt.grid()
plt.legend(loc= 'best')
plt.tight_layout()
plt.savefig(outDatabaseDir+'cell_average_elong_ROI_comparison.png')
```


![png](output_47_0.png)

