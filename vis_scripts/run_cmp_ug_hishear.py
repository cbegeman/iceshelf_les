# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np

sys.path.append('/Users/cbegeman/Software_files/my_python_code/')
base_dir = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'
run1 = 'test_ocean_amb_hishear_medres'
run2 = 'test_ocean_ug_hishear_medres'
#run2 = 'test_ocean_amb_hisheary_hiresz'
filedir1 = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'+run1+'/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/'
filedir2 = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'+run2+'/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/'
run = ['dpdx','dpdy']

plotvar_t = ['umax','vmax','melt','w"u"0','w"v"0','u*','dt','E','E*','ol','k_offset_mcph']
plotvar_pr = ['pt','sa','prho','rho_ocean','w',
              'w*u*','w"u"','w*v*','w"v"','vel_var_ratio',
              'e*','e','km','kh','l','wpt','wsa','w"pt"','w*pt*','w"sa"','w*sa*','wu','wv']
plotvar = ['pt','sa','rho_ocean','v','w','e']#'u',
plotvar2 = ['melt*_xy','shf*_xy','sasws*_xy','u*_xy','ol*_xy',
            'pt1*_xy','pt_io*_xy','thermal_driving','sa1*_xy','sa_io*_xy',
            'haline_driving','usws*_xy','vsws*_xy']

#tplot = np.arange(50,60,5)
tplot = [10.] 
#palm.plot_tseries_zlevel([filedir2],[run2],['sa','wu','U','wv'],zeval = [1.],tlim = [1.,20.])
#palm.plot_tseries_zlevel([filedir2],[run2],['sa'],zeval = [1.],tlim = [1.,20.],norm='mean')
#palm.plot_tseries([filedir1,filedir2], run, plotvar_t)
#palm.plot_tseries([filedir1,filedir2], [run1,run2], ['melt'])
#palm.plot_pr([filedir1,filedir2],run,plotvar = plotvar_pr,teval = tplot,norm=['' for i in plotvar_pr],zlim=[-500,0])
palm.plot_pr_varcmp([filedir1,filedir2],run,varcmp = ['velocity','variance','momflux_u','momflux_v','heatflux_z','saltflux_z'],teval = [tplot[0] for i in run])
#    palm.plot_pr_varcmp([filedir1,filedir2],run,varcmp = ['velocity','variance','momflux_u','momflux_v','heatflux_z','saltflux_z'],teval = [j])
#    palm.plot_pr_varcmp([filedir1,filedir2],run,varcmp = ['variance'],teval = [j])
#    for idx,i in enumerate([filedir2]):
#        palm.plot_uv_vector(i,run[idx],zeval = np.arange(0,-2*512/3,-10), teval = [j])
#for idx,i in enumerate([filedir1,filedir2]):
#    palm.plot_pr([i],[run[idx]],plotvar = plotvar_pr,teval = np.arange(3.,15.,1.))
#palm.plot_pr([filedir2],[run2],plotvar = ['N2'],teval = np.arange(6.,21.,1.),zlim=[-100,0])
    
#   palm.plot_3d_slice(diri,run[i],plotvar,teval = tplot, slice_dim = 'y')

#palm.plot_3d_slice(filedir1,run1,['w'],teval = tplot, slice_dim = 'z', zeval=[-300.,-100.,-10])
#palm.plot_2d_xy(filedir1,runname1,plotvar2,teval = np.arange(20,21,1))
#palm.plot_pr_TS(filedir1,run1,zeval = np.arange(0,-2*512/3,-10), teval = tplot)
