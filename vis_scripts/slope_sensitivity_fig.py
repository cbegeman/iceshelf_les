# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize2
from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

palm.plot_tseries_cross(diri, run, runlabel, 
                        plotvar = ['sin_alpha','melt'], 
                        data_type = ['parameter','ts'],
                        tav=tav_ts, teval=[tprofile,tprofile],
                        plot_legend=True, legtitle = legtitle,
                        figsize = figsize2, plot_cycles=True,
                        col = colorVal, overwrite=True)
palm.plot_tseries_cross(diri, run, runlabel, 
                        plotvar = ['u*','gamma_T_2m'], 
                        data_type = ['ts','pr'],
                        tav=tav_ts, teval=[tprofile,tprofile],
                        plot_legend=False, legtitle = legtitle,
                        figsize = figsize2, plot_cycles=True,
                        col = colorVal, overwrite=True)
