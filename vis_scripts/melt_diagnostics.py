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

palm.plot_tseries(diri, run, ['thermal_driving_i'], runlabel = runlabel, legtitle = legtitle, col=colorVal, figsize = figsize3,tlim = [tmin,tmax])
palm.plot_tseries_cross(diri, run, runlabel, 
                        plotvar = ['thermal_driving_i','melt'], 
                        data_type = ['ts','ts'],
                        tav=tav_ts, teval=[tprofile,tprofile],
                        plot_legend=False, legtitle = legtitle,
                        figsize = figsize2, plot_cycles=True, 
                        col = colorVal, overwrite=True)
palm.plot_tseries_cross(diri, run, runlabel, 
                        plotvar = ['u*','melt'], 
                        data_type = ['ts','ts'],
                        tav=tav_ts, teval=[tprofile,tprofile],
                        plot_legend=False, legtitle = legtitle,
                        figsize = figsize2, plot_cycles=True, 
                        col = colorVal, overwrite=True)
palm.plot_tseries_cross(diri, run, runlabel, 
                        plotvar = ['gamma_T','melt'], 
                        data_type = ['ts','ts'],
                        tav=tav_ts, teval=[tprofile,tprofile],
                        plot_legend=False, legtitle = legtitle,
                        figsize = figsize2, plot_cycles=True, 
                        col = colorVal, overwrite=True)
