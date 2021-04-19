#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 07:52:34 2019

@author: cbegeman
"""
import math,os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as tick
import my_array_mod as ar
import my_string_mod as strop
import convert_NCAR_mod as ncar
import extract_var_palm as pv 

filedir = 
uvar = pv.compute_variance(filedir,'u')
