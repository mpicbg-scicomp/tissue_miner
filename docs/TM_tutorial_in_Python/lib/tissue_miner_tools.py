import sqlite3 as lite
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import pandas.io.sql as psql
import pandas.io.parsers as pp

from sets import Set
import os.path
import time


def subset_dataframe(df, columns, conditions):
    """
    Subset a DataFrame object on columns in columns list
    by values is conditions list
    """ 
    t= df.copy()
    for column, condition in zip(columns, conditions):
        t= t[t[column] == condition]
    return t                        

def fill_zeros(s,n):
    """
    Add zeros to a string s until
    it reaches length n.
    """
    while len(s) < n:
        s= ''.join(('0',s))
    return s

def smooth_data(x, NSmooth=10, mode= 'valid'):
    """
    Returns the original data smoothed over a
    window of lenght NSmooth.
    """
    return np.convolve(x, 1.*np.ones(NSmooth)/NSmooth, mode= mode)


def make_directory(path):
    """
    Create all directories in a proposed path
    - unlike simple os.makedirs(path) for which all
    higher level directories have to exist.
    """
    path_list= path.rstrip().rstrip('/').split('/')
    making_path= ''
    if path_list[0] == '':
        making_path+= '/'
        path_list.pop(0)
    for segment in path_list:
        making_path+= segment
        if not os.path.exists(path):
            os.makedirs(path)
        making_path+='/'

def cells_total_area(rc):
    """
    Calculates sum of cell areas in each frame. 
    Intednded for use on 'cells' or subsets of 'cells' tables.
    """    
    return rc[['frame', 'area']].groupby('frame').agg(np.sum).reset_index().sort('frame')['area'].values

def cells_average_area(rc):
    """
    Calculates average of cell areas in each frame. 
    Intednded for use on 'cells' or subsets of 'cells' tables.
    """
    return rc[['frame', 'area']].groupby('frame').agg(np.mean).reset_index().sort('frame')['area'].values

def cells_average_elongation(rc):
    """
    Calculates average _xx and _xy components of elongation nematic in each frame.
    Intedend for use on 'cells' or subsets of 'cells' tables.
    """
    Q_xx= rc[['frame', 'elong_xx']].groupby('frame').agg(np.mean).reset_index().sort('frame')['elong_xx'].values
    Q_xy= rc[['frame', 'elong_xy']].groupby('frame').agg(np.mean).reset_index().sort('frame')['elong_xy'].values
    return Q_xx, Q_xy

def cells_number(rc):
    """
    Calculates number of cells in each frame.
    Intedend for use on 'cells' or subsets of 'cells' tables.
    """    
    return rc[['frame', 'cell_id']].groupby('frame').agg(len).reset_index().sort('frame')['cell_id'].values


