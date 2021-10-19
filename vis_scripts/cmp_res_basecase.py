# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
import numpy as np
import matplotlib as pltlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10
from run_table_mod import value_from_namelist
from var_param_palm import K0
import cmocean as cmo
import run_table_mod as table
from var_param_palm import *
from plot_param_palm import *
import plot_palm_mod as palm
from gsw import t_freezing, pt0_from_t

pltlib.rc_file('rcparams.txt', use_default_template=True)
sens_var = 'pt'

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')
base_dir = ['/turquoise/usr/projects/climate/cbegeman/palm/jobs/',
            '/lustre/scratch3/turquoise/cbegeman/palm/jobs/']
filedir = ['test_ocean_melt_batch_dpdy_021','test_ocean_melt_batch_dpdy_021_hires']
rundir='/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/'
diri = [base_dir[i]+filedir[i]+rundir for i in range(2)]
legtitle=r'Resolution (m)'
linestyle = ['-']

runlabel = [r'\Delta x,y=0.5, \Delta z=0.25',
            r'\Delta x,y=0.25, \Delta z=0.125']
run = ['dz1, dz0.5']
tend = np.zeros((len(filedir),))
for idx,i in enumerate(diri):
    runfile = i + 'RUN_CONTROL'
    tend[idx] = palm.end_time(runfile)/3600.

tplot = np.min(tend)
tcross = np.min(tend)

tmin = 2.
if tplot <2:
   tmin = 0.
tunits = 'hr'
tmax = 50
tav_pr = 13
tav_ts = 13.
tprofile = tplot-tav_pr/2

xmid = 64.
ymid = 64.
zlim = [-5,1.]
zmax = -2*64/3
z_cfl = -2

#plot_tseries  = False
plot_cross    = False
#plot_profiles = False
plot_slices   = False
plot_hovmoller = False

plot_tseries   = True
#plot_cross     = True
plot_profiles  = True
#plot_slices    = True
#plot_hovmoller = True

tidal=False
overwrite=False
