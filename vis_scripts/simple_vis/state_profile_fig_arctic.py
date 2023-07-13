# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize3
#from cmp_thermal_driving_slope1 import *
from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

#palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [tplot[0]], ops=norm_pr,zlim=[zmax,0])
palm.plot_pr(diri, run, ['pt'],
             runlabel=runlabel,  
             teval = [tprofile,tprofile], 
             ops=['far_diff'], tav = tav_pr, 
             figsize = figsize3,legtitle=legtitle,
             #xscale_input = xscale_input, xscale_label = xscale_label, zscale = 'Ekman', 
             zlim=[zmax,0], col=colorVal)#, xlim = [0,20])
palm.plot_pr(diri, run, ['sa','velocity'],
             runlabel=runlabel,  
             teval = [tprofile,tprofile], 
             ops=['far_diff',''], tav = tav_pr, 
             figsize = figsize3,
             #xscale_input = xscale_input, xscale_label = xscale_label, zscale = 'Ekman', 
             plot_legend = False,
             zlim=[zmax,0], col=colorVal)#, xlim = [0,20])
#palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [min(tplot),max(tplot)],tall=True, ops=norm_pr,zlim=[zmax,0])
