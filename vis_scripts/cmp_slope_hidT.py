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
from math import log10,pi
from run_table_mod import value_from_namelist
from var_param_palm import K0
import cmocean as cmo
import run_table_mod as table
from var_param_palm import *
from plot_param_palm import *
import plot_palm_mod as palm

pltlib.rc_file('rcparams.txt', use_default_template=True)

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')
#base_dir = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'
base_dir = '/turquoise/usr/projects/climate/cbegeman/palm/jobs/'
filedir = ['test_ocean_melt_batch_alpha_surface_0{:02d}'.format(i) for i in range(6,9)]
filedir.append('test_ocean_melt_batch_dpdy_021')
diri = [base_dir+i+'/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/' for i in filedir]
legtitle=r''

sens_var = 'alpha_surface'
linestyle = ['-']

run= ['' for i in filedir]
runlabel = ['' for i in filedir]
slope = np.zeros((len(filedir),))
for idx,i in enumerate(filedir):
    slope[idx] = value_from_namelist(base_dir+i,'alpha_surface')[0]
    runlabel[idx] = '{:2.1e}'.format(slope[idx])
    run[idx] = 'slope{:1.0e}'.format(slope[idx])
#xscale_input = np.sin(pi*slope/180.)
#xscale_input = np.power(np.sin(pi*slope/180.),0.5)
xscale_input = np.power(np.sin(pi*slope/180.),0.25)
#xscale_input = np.log(np.sin(pi*slope/180.))
print(xscale_input)
xscale_label = r'/ \sin\alpha^{1/4} \:('

cNorm  = colors.Normalize(vmin=np.min(np.log10(slope))-0.5,vmax=np.max(np.log10(slope))+0.5)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='cmo.ice_r')
colorVal = np.zeros((len(slope),4))
for idx,i in enumerate(slope):
    colorVal[idx,:] = scalarMap.to_rgba(np.log10(i))

a = np.array([[0,1]])
plt.figure(figsize=(6, 3))
img = plt.imshow(a, norm=cNorm, cmap='cmo.ice_r')
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
cbar = plt.colorbar(orientation="horizontal", cax=cax)
cbar.ax.tick_params(labelsize=16)
cbar_label = r'$log(\alpha) (^{\circ})$'
cbar.set_label(cbar_label,fontsize = 20)
plt.savefig('colorbar_slope.png')

tend = np.zeros((len(filedir),))
for idx,i in enumerate(diri):
    runfile = i + 'RUN_CONTROL'
    tend[idx] = palm.end_time(runfile)/3600.
print(runlabel)
print(tend)

tperiod = 12.
tmin = 2.
tunits = 'hr'
tav_pr = 13.
tav_ts = 12.
tmax = 52.
tplot = 50#np.min(tend)
#tplot = np.max(tend)
#tcross = np.min(tend)
print('tplot=',tplot)
tcross = tmin
tprofile = tplot-tav_pr/2#np.arange(39.,50.,1.)#[39.,50.]#[np.min(tend)]

ymid = 64.
xmid = 64.
zmax = -2*64/3
z_cfl = -2

plot_tseries  = False
plot_cross    = False
#plot_profiles = False
plot_slices   = False
#plot_tseries   = True
#plot_cross     = True
plot_profiles  = True
#plot_slices    = True

tidal=False
overwrite=False
#palm.plot_tseries(diri,run,['melt'],tlim=[tmin,tplot],
#                  runlabel=runlabel,legtitle=legtitle,
#                  col = colorVal, 
#                  linestyle=['None'],marker='.',
#                  overwrite=overwrite)
