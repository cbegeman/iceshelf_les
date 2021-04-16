"""
Created on Tue Jan 22 12:50:47 2019

@author: cbegeman
"""

import numpy as np

def last_index(self):
    return len(self)-1;

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    val = array[idx]
    return val,idx;
