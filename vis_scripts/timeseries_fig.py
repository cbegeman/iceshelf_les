# -*- coding: utf-8 -*-
"""
Timeseries figure

"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize3
#from cmp_thermal_driving_slope1 import *
from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

palm.plot_tseries(diri, run, ['E*'], runlabel = runlabel, legtitle = legtitle, col=colorVal, figsize = figsize3,tlim = [tmin,tmax], plot_legend=False, ylim=[4e-7,3e-5])
palm.plot_tseries(diri, run, ['melt'], runlabel = runlabel, legtitle = legtitle, col=colorVal, figsize = figsize3,tlim = [tmin,tmax], plot_legend=False, linestyle = ['None' for i in diri], marker = '.', ylim=[2e-3,7])
palm.plot_tseries(diri, run, ['u*'], runlabel = runlabel, legtitle = legtitle, col=colorVal, figsize = figsize3,tlim = [tmin,tmax], ylim=[0,6e-3])
