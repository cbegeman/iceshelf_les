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
#run1 = 'test_ocean_bsign_noslp_dSdz2_medres'
#run1 = 'test_ocean_bsign_noslp_lowres'
#run1 = 'test_ocean_bsign_modslp_warm_medres'
run1 = 'test_ocean_bsign_modslp_hiresz'
filedir1 = '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'+run1+'/RUN_ifort.grizzly_hdf5_mpirun_test_oceanml/'

plotvar_t = ['umax','vmax','melt','w"u"0','w"v"0','u*','dt','E','E*','ol','k_offset_mcph']
plotvar = ['pt','sa','u','v']
plotvar_pr = ['pt','sa','prho','w','Sw','S2','u2v2','u','v',
              'w*u*','w"u"','w*v*','w"v"','vel_var_ratio',
              'e*','e','km','kh','l','wpt','wsa','w"pt"','w*pt*','w"sa"','w*sa*','wu','wv']
plotvar_varcmp = ['velocity','variance','momflux_z','heatflux_z','saltflux_z','tke']
plotvar_flux = ['wv','wu','wpt','wsa']
norm_pr = ['' for i in plotvar_pr]
norm_pr[plotvar_pr.index('pt')] = 'far_diff'
norm_pr[plotvar_pr.index('sa')] = 'far_diff'
norm_3d = ['' for i in plotvar_pr]
norm_3d[plotvar_pr.index('sa')] = 'mean'
plotvar = ['pt','sa','u','v','w','e']
plotvar2 = ['melt*_xy','u*_xy','ol*_xy',
            'pt1*_xy','pt_io*_xy','thermal_driving','sa1*_xy','sa_io*_xy',
            'haline_driving']#,'usws*_xy','vsws*_xy']
plotvar_zlevel = ['pt','sa','u','v']
runname = ''
zmax = -50
#zmax = -500

#tplot = [98.]
tplot = [6.]
#tplot = np.arange(40,99,5)
#tplot = np.arange(5,15,5)

#palm.plot_tseries([filedir1], [runname], plotvar_t)
#palm.plot_tseries_melt_us(filedir1, runname)
#palm.plot_tseries_zlevel([filedir1], [runname], plotvar_zlevel,norm=['' for i in plotvar_zlevel], zeval=[-10.], tlim=[tplot[0],tplot[-1]])

#palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [tplot[0]], ops=norm_pr,zlim=[zmax,0])
#palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [min(tplot),max(tplot)],tall=True, ops=norm_pr,zlim=[zmax,0])
#palm.plot_pr_varcmp([filedir1],[runname],plotvar_varcmp,teval = tplot,zlim=[zmax,0])

#palm.plot_uv_vector(filedir1,run1,zval = [-2*512/3,0], teval = tplot)
#palm.plot_pr_TS(filedir1,run1,zval = [-2*512/3,0], teval = tplot)
#palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot,slice_dim = 'y',ops=['' for i in plotvar],xval=[125.,125.],zval=[zmax,0])
palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot,slice_dim = 'z',ops=['' for i in plotvar],zval=[-1,-1])
#palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot,slice_dim = 'y',ops=['' for i in plotvar],yval=[125.,125.],zval=[zmax,0])
#palm.plot_3d_slice(filedir1,runname,plotvar = ['u"'],teval = tplot,ops=[''],slice_dim='y')
#for j in [-1.]:
#    palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot, slice_dim = 'z', zeval=[j,j],norm=norm_3d)
#palm.plot_2d_xy(filedir1,runname,['melt*_xy'],teval = [10.])
