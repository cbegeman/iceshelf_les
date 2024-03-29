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

palm.plot_pr(diri, run, ['vel_var_ratio'],#,'km_eff'
             runlabel=runlabel,  
#             teval = [40,40], 
             teval = [tprofile,tprofile], 
             tav = tav_pr,
             figsize = figsize2, 
             legtitle=legtitle,
             xlim = [0,1], zlim=[-20,0], col=colorVal
            )
#palm.plot_pr(diri, run, ['wsa'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr,
#             figsize = figsize2, 
#             legtitle=legtitle,
#             zlim=[zmax,0], col=colorVal
#            )
#palm.plot_pr(diri, run, ['km_eff'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr,
#             figsize = figsize2, 
#             legtitle=legtitle,
#             zlim=[-24,0], col=colorVal
#            )
#palm.plot_pr(diri, run, ['k_all'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr, 
#             figsize = figsize2, 
#             legtitle=legtitle,
#             plot_legend=False,
#             zlim=[zmax,0], col=colorVal
#            )
