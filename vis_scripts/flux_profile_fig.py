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

#palm.plot_pr(diri, run, ['wpt','wsa'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr,
#             figsize = figsize2, 
#             legtitle=legtitle,
#             zlim=[zmax,0], col=colorVal
#            )
#palm.plot_pr(diri, run, ['momflux_z'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr, 
#             figsize = figsize2, 
#             legtitle=legtitle,
#             plot_legend=False,
#             zlim=[zmax,0], col=colorVal
#            )
#palm.plot_pr(diri, run, ['wpt','momflux_z'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr, 
#             xscale_input = xscale_input, 
#             xscale_label = xscale_label, 
#             zscale = 'Ekman', 
#             figsize = figsize2, 
#             plot_legend=False,
#             #zlim=[-25.,0.], 
#             col=colorVal
#            )
#palm.plot_pr(diri, run, ['heatflux_z'],#'saltflux_z','momflux_u','momflux_v'],
#             runlabel=runlabel,  
#             teval = [tprofile,tprofile], 
#             tav = tav_pr,
#             figsize = figsize2, 
#             legtitle=legtitle,
#             plot_legend=False,
#             zlim=[zmax,0], col=colorVal
#            )
palm.plot_pr(diri, run, ['heatflux_z'],#'saltflux_z','momflux_u','momflux_v'],
             runlabel=runlabel,  
             teval = [tmin+tav_pr/2,tmin+tav_pr/2], 
             tav = tav_pr,
             figsize = figsize2, 
             legtitle=legtitle,
             plot_legend=False,
             zlim=[zmax,0], col=colorVal
            )
