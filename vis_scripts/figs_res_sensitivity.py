# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from cmp_res_basecase import *

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')

plotvar_t = ['melt']
plotvar_pr = ['sa','velocity','heatflux_z']
norm_pr = ['' for i in plotvar_pr]
if plot_tseries:
    palm.plot_tseries(diri, run, plotvar_t, 
                      runlabel=runlabel, tlim=[tmin,tplot], 
                      linestyle=['None'], linewidth=0.,
                      marker='.',
                      figsize=figsize2)
if plot_profiles:
    palm.plot_pr(diri, run, ['pt'],
                 runlabel=runlabel,  
                 teval = [tprofile,tprofile], 
                 ops=norm_pr, tav = tav_pr, 
                 figsize=figsize3, zlim=[zmax,0], legtitle=legtitle)
    palm.plot_pr(diri, run, plotvar_pr,
                 runlabel=runlabel,  
                 teval = [tprofile,tprofile], 
                 ops=norm_pr, tav = tav_pr, 
                 figsize=figsize3, zlim=[zmax,0], plot_legend=False)
