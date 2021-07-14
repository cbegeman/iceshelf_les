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
from cmp_thermal_driving_slope1 import *
#from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

palm.plot_hovmoller(diri,run,['e*'],
                    tlim = [tmin,tmax], runlabel=runlabel,zlim=[zmax,0],overwrite=True, clim=[-7.5,-4.25], 
                    figsize = figsize2_box
#                    #contour_var = 'vel_var_ratio', contour_val = 0.1
                    ,contour_var = 'diss', contour_val = 1e-9
                   )
#palm.plot_hovmoller(diri,run,['u*2','v*2'],
#                    tlim = [tmin,tmax], runlabel=runlabel,zlim=[zmax,0],overwrite=True, 
#                    figsize = figsize2_box
#                    #contour_var = 'vel_var_ratio', contour_val = 0.1
#                   )
#palm.plot_hovmoller(diri,run,['diss'],
#                    tlim = [tmin,tmax], runlabel=runlabel,zlim=[zmax,0],overwrite=True, 
#                    figsize = figsize2_box
#                    ,contour_var = 'diss', contour_val = 1e-9, clim=[1e-9,1e-7]
#                   )
#palm.plot_hovmoller(diri,run,['vel_var_ratio'],tlim = [tmin,tmax], runlabel=runlabel,zlim=[zmax,0],overwrite=True, clim=[0.1,2])
