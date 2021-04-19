# -*- coding: utf-8 -*-
"""
"""

import sys
import numpy as np
import matplotlib as pltlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10
import run_table_mod as table
from run_table_mod import value_from_namelist
import cmocean as cmo
import run_table_mod as table
from var_param_palm import *
from plot_param_palm import *
import plot_palm_mod as palm

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')
base_dir = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'
filedir = ['test_ocean_amd_mcphee_hires_slope1e-1_dpdx-3e-2_dpt1e0' 
          ,'test_ocean_amd_mcphee_hires_rot_slope1e-1_dpdx3e-2_dpt1e0' 
          ,'test_ocean_amd_mcphee_hires_slope1e-1_dpdy-3e-2_dpt1e0_new' 
          ,'test_ocean_amd_mcphee_hires_slope1e-1_dpdx--3e-2_dpt1e0' 
          ,'test_ocean_amd_mcphee_hires_slope1e-1_dpdy--3e-2_dpt1e0' 
          ,'test_ocean_amd_mcphee_hires_rot_slope1e-1_dpdy3e-2_dpt1e0' 
          ,'test_ocean_amd_mcphee_hires_rot_slope1e-1_dpdx-3e-2_dpt1e0' 
          ,'test_ocean_amd_mcphee_hires_rot_slope1e-1_dpdy-3e-2_dpt1e0'] 
linestyle = ['-','--','-','-','-','--','--','--']
diri = ['/lustre/scratch3/turquoise/cbegeman/palm/jobs/'+i+'/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/' for i in filedir]
#run = ['all_dpdx_3e-2','dpdy_3e-2','dpdx_-3e-2','dpdy_-3e-2']
legtitle=r''
#runlabel = ['dpdx = 3e-2','dpdy = 3e-2','dpdx = -3e-2','dpdy = -3e-2']

sens_var = ''
[varname,runname,runtable] = table.load_vals('run_table.txt')
run = ['' for i in filedir]
runlabel = ['' for i in filedir]
dpdx = np.zeros((len(filedir),))
dpdy = np.zeros((len(filedir),))
for idx,i in enumerate(filedir):
    dpdx[idx] = runtable[runname.index(i),varname.index('dpdxy')]
    dpdy[idx] = runtable[runname.index(i),varname.index('dpdxy_y')]
    runlabel[idx] = 'dpdx = [{:2.2f},{:2.2f}]'.format(dpdx[idx],dpdy[idx])
    if dpdx[idx] != 0:
       run[idx] += 'dpdx_{:2.1e}'.format(dpdx[idx])
    if dpdy[idx] != 0:
       run[idx] += 'dpdy_{:2.1e}'.format(dpdy[idx])

run= ['' for i in filedir]
runlabel = ['' for i in filedir]
dpdx = np.zeros((len(filedir),))
dpdy = np.zeros((len(filedir),))
for idx,i in enumerate(filedir):
    varvalue = value_from_namelist('/lustre/scratch3/turquoise/cbegeman/palm/jobs/'+i,'dpdxy')
    dpdx[idx] = varvalue[0]
    dpdy[idx] = varvalue[1]
    runlabel[idx] = 'dpdx = [{:2.2f},{:2.2f}]'.format(dpdx[idx],dpdy[idx])
    if dpdx[idx] != 0:
       run[idx] += 'dpdx_{:2.1e}'.format(dpdx[idx])
    if dpdy[idx] != 0:
       run[idx] += 'dpdy_{:2.1e}'.format(dpdy[idx])

tend = np.zeros((len(filedir),))
for idx,i in enumerate(diri):
    runfile = i + 'RUN_CONTROL'
    tend[idx] = palm.end_time(runfile)/3600.

tplot = np.max(tend)
tcross = np.min(tend)
tprofile = [np.min(tend)]
tmin = 5.
tunits = 'hr'
tav_pr = 1.
tav_ts = 2.
colorVal = col
colorVal.append(col)
print(runlabel)
print(tend)

xmid = 64.
ymid = 64.
zmax = -2*32/3
#zmax = -10.
dz = 0.25
z_cfl = -2

overwrite=False
tidal=False
