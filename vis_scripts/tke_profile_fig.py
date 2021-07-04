# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
#from cmp_thermal_driving_slope1 import *
from cmp_slope_hidT import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

plotvar_pr = ['Fshear','Ftrans','Fbuoy_uw']
norm_pr = ['' for i in plotvar_pr]

palm.plot_pr(diri, run, ['e*'],
             runlabel=runlabel,  
             teval = [tprofile,tprofile], 
             ops=norm_pr, tav = tav_pr, legtitle=legtitle, 
             figsize = figsize2, zlim=[zmax,0], col=colorVal)#, xlim = [0,20])
#palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [tplot[0]], ops=norm_pr,zlim=[zmax,0])
#palm.plot_pr(diri, run, plotvar_pr,
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             ops=norm_pr, tav = tav_pr, 
#             plot_legend=False,
#             #xscale_input = xscale_input, xscale_label = xscale_label, zscale = 'Ekman', 
#             figsize = figsize2, zlim=[zmax,0], col=colorVal)#, xlim = [0,20])
#palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [min(tplot),max(tplot)],tall=True, ops=norm_pr,zlim=[zmax,0])
