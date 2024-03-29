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
#base_dir = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'
base_dir = '/turquoise/usr/projects/climate/cbegeman/palm/jobs/'
filedir = ['test_ocean_melt_batch_dpdy_021']
for i in range(18,21):
    filedir.append('test_ocean_melt_batch_dT_0{:02d}'.format(i))
diri = [base_dir+i+'/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/' for i in filedir]
#legtitle=r'$T_{\infty} - T_f(S_{\infty})\:(^{\circ}C)$'
legtitle=r'$\Delta \theta_{\infty} (^{\circ}C)$'
xscale_input = np.ones((len(diri)))
xscale_label = r'/\Delta \theta_{\infty} \: (m\:s^{-1})$'
#print(xscale_input)
linestyle = ['-']

runlabel = ['' for i in filedir]
run = ['' for i in filedir]
pt_surf = np.zeros((len(filedir),))
tf = t_freezing(35,800,0)
pt_freezing = pt0_from_t(35,tf,800)
for idx,i in enumerate(filedir):
    pt_surf[idx] = value_from_namelist(
                   base_dir+i,
                   'pt_surface')[0] - K0
    thermal_driving = pt_surf[idx] - pt_freezing
    #runlabel[idx] = '{:2.1f}'.format(pt_surf[idx])
    runlabel[idx] =  '{:2.1f}'.format(thermal_driving)
    xscale_input[idx] = thermal_driving
    #runlabel[idx] = r'$T_{\infty} - T_f(S_{\infty} = $' + '{:2.1f}'.format(thermal_driving)
    run[idx] = r'dT'+runlabel[idx]

cNorm  = colors.Normalize(vmin=np.min(pt_surf)-0.05,vmax=np.max(pt_surf)+0.1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='cmo.thermal')
colorVal = np.zeros((len(pt_surf),4))
for idx,i in enumerate(pt_surf):
    colorVal[idx,:] = scalarMap.to_rgba(i)

#a = np.array([[0,1]])
#plt.figure(figsize=(6, 3))
#img = plt.imshow(a, norm=cNorm, cmap='cmo.thermal')
#plt.gca().set_visible(False)
#cax = plt.axes([0.1, 0.2, 0.8, 0.6])
#cbar = plt.colorbar(orientation="horizontal", cax=cax)
#cbar.ax.tick_params(labelsize=16)
##cbar_label = r'Far-field temperature $^{\circ}C$'
#cbar_label = r'$\theta ^{\circ}C$'
#cbar.set_label(cbar_label,fontsize = 20)
#plt.savefig('colorbar_thermal_driving.png')

tend = np.zeros((len(filedir),))
for idx,i in enumerate(diri):
    runfile = i + 'RUN_CONTROL'
    tend[idx] = palm.end_time(runfile)/3600.

#print(filedir)
#print(runlabel)
#print(colorVal)
#print(tend)

tplot = 50.
#tplot = np.min(tend)
#tplot = np.max(tend)
tcross = np.min(tend)
tperiod = 12.

tmin = 2.
if tplot <2:
   tmin = 0.
tunits = 'hr'
tmax = 50
#tav_pr = 1
tav_pr = 13
tav_ts = 13.
tprofile = tplot-tav_pr/2#[40.,48.]

xmid = 64.
ymid = 64.
zlim = [-5,1.]
zmax = -2*64/3
z_cfl = -2

plot_tseries  = False
plot_cross    = False
plot_profiles = False
#plot_slices   = False
plot_hovmoller = False
#plot_tseries   = True
#plot_cross     = True
#plot_profiles  = True
plot_slices    = True
#plot_hovmoller = True

tidal=False
overwrite=False
#palm.plot_tseries(diri,run,['melt'],tlim=[tmin,tplot],
#                  runlabel=runlabel,legtitle=legtitle,
#                  col = colorVal, 
#                  linestyle=['None'],marker='.',
#                  overwrite=overwrite)
