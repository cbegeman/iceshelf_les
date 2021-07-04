# -*- coding: utf-8 -*-
"""
Timeseries figure

"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize2
#from cmp_thermal_driving_slope1 import *
from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

palm.plot_tseries(diri, run, ['w*2'], runlabel = runlabel, legtitle = legtitle, col=colorVal, figsize = figsize2,tlim = [tmin,tmax], linestyle = ['None' for i in diri], marker = '.', ylim=[2e-6,3e-5])
