import matplotlib as mpl
import sqlite3 as lite
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp
import rpy2.robjects as robjects
import rpy2.robjects as ro

import matplotlib.cm as cm
from matplotlib import collections as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sets import Set
import os.path

import time

import tissue_miner_tools as tml

class Movie:
    """
    Movie class is intended to contain all the data of a single time-lapse movie.

    It also contains methods for loading the data from SQL database and extra .RData
    tables provided by TissueMiner workflow.

    Note that it will create pickle (.pkl) format of tables provided in RData files
    to significantly reduce data loading time. Currently this is done trough intermediate
    step creating .csv file, which seems to be faster than using rpy2 library to
    directly read RData files into python for large files.
    Each time SQL database and/or RData files are updated piclke files will be updated.
    """
    def __init__(self, name, path, ROI_path= '', ROI_filename= ''):
        self.name= name
        self.DB_path= path
        self.con= lite.connect(self.DB_path+name+'/'+name+'.sqlite')
        self.loaded= Set()
        self.DB_table= {}
        self.intermediate_state= {}
        self.triangle_state= {}
        self.pupalWing_loaded= False
        try:
            print('Loading table frames from database ' + self.name + '...')             
            self.DB_table['frames']= psql.read_sql("SELECT * FROM frames;", self.con)
            self.loaded.add('frames')
            self.load_DB_table('cells')
            self.frames= self.DB_table['cells']['frame'].unique()
            self.time= self.DB_table['frames']['time_sec'].values/3600.
            self.dt= self.time[1:] - self.time[:-1]
            self.NFrames= len(self.frames)
        except:
            print('Table frames not available in ' + self.name + ' database.')
        try:
            pupalWings= pp.read_csv(self.DB_path+'PupalWingMovies.csv', sep=',')
            self.pupalWing_loaded= True
        except:
            print('PupalWingMovies.csv not found in: '+self.DB_path)
        try:
            self.time_shift= float(np.array(pupalWings[pupalWings['name']==name]['time_shift_sec'])[0])/3600.
            if np.isnan(self.time_shift):
                print('While loading ' + name + ' time shift is NAN')
            else:    
                self.time+= 15. + self.time_shift
        except:
            print('While loading '+ name+'. No time shift provided!.')  
        try:
            self.load_roiBT(path= ROI_path, filename= ROI_filename)            
        except:
            print('ROI file is not available!')
            
    def load_DB_table(self, table_name):
        """
        Loads a table table_name from the SQL database associated with the class
        instance on initialization.
        If pickle format of the table exists and is newer than the database file,
        it will load the table from pickle file instead.

        Also checks whether the table is already loaded.
        """
        if not table_name in self.loaded:
            file_pickle= self.DB_path + self.name + '/' + table_name + '.pkl'
            file_sql= self.DB_path + self.name + '/' + self.name + '.sqlite'
            if ((not os.path.isfile(file_pickle)) or (os.path.getmtime(file_sql) > os.path.getmtime(file_pickle))):
                print('Loading table ' + table_name + ' from database ' + self.name + '...')
                if table_name in ('vertices', 'bonds', 'frames') :
                    self.DB_table[table_name]= psql.read_sql('SELECT * FROM ' + table_name + ';', self.con)
                else:
                    self.DB_table[table_name]= psql.read_sql('SELECT * FROM ' + table_name + ' WHERE cell_id > 10000;', self.con)

                print('Writing table ' + table_name + ' to pickle file ' + self.name)
                self.DB_table[table_name].to_pickle(file_pickle)
            else:
                print('Loading table ' + table_name + ' from pickle file...') 
                self.DB_table[table_name]= pd.read_pickle(file_pickle)
            if table_name == 'cells':
                self.DB_table[table_name]= self.DB_table[table_name][['frame', 'cell_id', 'center_x', 'center_y', 'area', 'elong_xx', 'elong_xy']]                
            self.loaded.add(table_name)
            
    def database_tables(self):
        """
        Returns a list of tables of the SQL database as DataFrame.
        """
        return psql.read_sql("SELECT name FROM sqlite_master WHERE type='table';", self.con)
    
    def RData_to_pickle(self, table, file_RData, file_pickle):
        """
        Writes a table from RData file to a pickle format if RData file is newer
        that the pickle file, or pickle file does not exist.
        """
        file_csv= file_RData[:-6] + '.csv'
        if ((not os.path.isfile(file_pickle)) or
            (os.path.getmtime(file_RData) > os.path.getmtime(file_pickle))):
            print('Converting \n'+ file_RData + ' to \n' + file_pickle)
            print('Creating temporary .csv file...')            
            ro.r('load("'+file_RData+'")')
            ro.r('write.csv('+table+', "' + file_csv + '")')
            print('.csv created, writing to pickle')
            print('Reading temporary .csv file...')
            temp_csv= pp.read_csv(file_csv)
            print('Writing to pickle...')
            temp_csv.to_pickle(file_pickle)
        
    def load_roiBT(self, path, filename):
        """
        Loads ROI table.
        """
        if not 'roiBT' in self.loaded:
            file_pickle= self.DB_path + self.name + '/' + path + filename + '.pkl'
            file_RData= self.DB_path + self.name + '/' + path + filename + '.RData'
            print('Loading roiBT ...')
            self.RData_to_pickle(filename, file_RData, file_pickle)
            self.roiBT= pd.read_pickle(file_pickle)
            self.loaded.add('roiBT')
            self.regions= self.roiBT['roi'].unique()

    def load_cellshapes(self):
        """
        Loads cellshapes table.
        """
        if not 'cellshapes' in self.loaded:
            file_pickle= self.DB_path + self.name + '/cellshapes.pkl'
            file_RData= self.DB_path + self.name + '/cellshapes.RData'
            print('Loading cellshapes of ' + self.name + '...')
            self.RData_to_pickle('cellshapes', file_RData, file_pickle)
            self.cellshapes= pd.read_pickle(file_pickle)
            self.loaded.add('cellshapes')
            del self.cellshapes[self.cellshapes.columns[0]]            

    def show_image(self, frame):
        """
        Shows original (rotated) image of the given frame.
        """
        im_path= self.DB_path + self.name + '/Segmentation/' + self.name + '_' + tml.fill_zeros(str(frame), 3) + '/original_trafo.png'
        im= plt.imread(im_path)
        plt.imshow(im)
                                    
    def region_cells(self, region):
        """
        Returns cells DataFrame with cells belonging to the region ROI.
        """
        self.load_DB_table('cells')
        if region not in self.regions:
            raise Exception('Region '+region+' is not defined in this movie!')
        else:
            return self.DB_table['cells'][self.DB_table['cells']['cell_id'].isin(self.roiBT[self.roiBT['roi']==region]['cell_id'])]
                
    def region_deform_tensor(self, roi_name):
        """
        Loads precalculated deformation data.
        """
        if 'avgDeformTensorsWide.tsv' in os.listdir(self.DB_path+self.name+'/shear_contrib/'+roi_name):
            df_DB_shear= pp.read_csv(self.DB_path+self.name+'/shear_contrib/'+roi_name+'/avgDeformTensorsWide.tsv', sep='\t')
        return df_DB_shear

    def subset_table_by_region(self, df, region= 'blade', on_column= 'cell_id'):
        """
        Returns subset of DataFrame table with the rows corresponding to
        'cell_id' values in a given region.
        """
        rc= self.region_cells(region)['cell_id'].unique()
        return df[df[on_column].isin(rc)]

    def cell_number(self, region= 'whole_tissue'):
        """
        Calculates number of cells in a region for each frame.
        """
        region_cells= self.region_cells(region)
        nr_cells= region_cells[['frame', 'cell_id']].groupby('frame').agg(len).reset_index().rename(columns= {'cell_id': 'nr_cells'})
        return nr_cells

    def plot_frame_cells(self, frame, coll_df, color_column, c_min= 0., c_max= 1., n_ticks= 5, figsize= (6, 10), polygon_lw= .1, color_map= cm.afmhot, title= ''):
        """
        Plots a collection of polygons provided in coll_df DataFrame in 'plot_vertices' column.
        Color is assigned based on values in color_column column of the coll_df DataFrame.
        c_min and c_max control the range of the colormap.
        Colormap can be provided by user and is set to afmhot by default.
        """
        plt.figure(figsize= figsize)
        plt.title(title, fontsize= 25)
        self.show_image(frame)
        plt.gca().autoscale_view()
        plt.gca().set_aspect('equal')
        colors= color_map((coll_df[color_column].values-c_min)/(c_max - c_min)) 
        coll= mc.PolyCollection(coll_df['plot_vertices'].values, lw= polygon_lw)
        coll.set(facecolors= colors)
        plt.gca().add_collection(coll)
        plt.xlim(0, 900)
        plt.ylim(300, 1800)
        plt.gca().invert_yaxis()
        plt.axis('off')
        divider= make_axes_locatable(plt.gca())
        cax= divider.append_axes('right', size= '5%', pad= 0.05)
        mm= cm.ScalarMappable(cmap= color_map)
        mm.set_array(colors)
        cbar= plt.colorbar(mm, cax= cax, cmap= color_map, ticks= np.linspace(0, 1, n_ticks + 1))
        cbar.ax.set_yticklabels(np.linspace(c_min, c_max, n_ticks + 1))
