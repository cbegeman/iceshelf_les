# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize2_wide,figsize2_box
#from cmp_thermal_driving_slope1 import *
from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

for i,filedir in enumerate([diri[0],diri[-1]]):
     palm.plot_3d_slice(filedir,run[i],['w'],
                        teval = [40],
                        slice_dim = 'y',clim=[-2.44,-2.39],
                        zval=[zmax,0],yval=[64,64],
                        figsize=figsize2_wide)
     palm.plot_3d_slice(filedir,run[i],['sa'],
                        teval = [40],
                        slice_dim = 'y',clim=[34.96,35],
                        zval=[zmax,0],yval=[64,64],
                        figsize=figsize2_wide)
     palm.plot_3d_slice(filedir,run[i],['u','v'],
                        teval = [40],
                        slice_dim = 'z',clim=[-2e-2,2e-2], 
                        zval=[-1,-1],
                        figsize=figsize2_box)
#        for z in [-1]:#,-5,-10,-20]:
#        #for z in [-1,-5,-10,-20]:
#            palm.plot_3d_slice(filedir,run[i],plotvar,teval = np.arange(40,49,1),slice_dim = 'z',zval=[z,z])
#palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot,slice_dim = 'y',ops=['' for i in plotvar],yval=[125.,125.],zval=[zmax,0])
#palm.plot_3d_slice(filedir1,runname,plotvar = ['u"'],teval = tplot,ops=[''],slice_dim='y')
#palm.plot_2d_xy(filedir1,runname,['melt*_xy'],teval = [10.])
