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

plotvar_t = ['umax','vmax','melt','w"u"0','w"v"0','u*','dt','E','E*','ol','k_offset_mcph']
plotvar = ['pt','sa']#,'u','v']
plotvar_pr = [#'pt','sa','prho','w','Sw','u2v2','u','v',
              'N2','S2'#,'Rf'
              #'k_all','km','kh','km_eff',
              #,'velocity','momflux_z','w*u*','w"u"','w*v*','w"v"','wu','wv'
              #,'vel_var_ratio'
              ,'Fshear','Ftrans','Fbuoy_uw'#,'w*p*:dz','w*u*u*:dz','Fbuoy','Fbuoy_u','Fbuoy_w','tke_all'
              ,'e*','e','l'#,'wpt','wsa','w"pt"','w*pt*','w"sa"','w*sa*',
              ]
plotvar_varcmp = ['velocity','variance','momflux_z','heatflux_z','saltflux_z','tke']
plotvar_flux = ['wv','wu','wpt','wsa']
norm_pr = ['' for i in plotvar_pr]
#norm_pr[plotvar_pr.index('pt')] = 'far_diff'
#norm_pr[plotvar_pr.index('sa')] = 'far_diff'
#norm_3d = ['' for i in plotvar_pr]
#norm_3d[plotvar_pr.index('sa')] = 'mean'
plotvar2 = ['melt*_xy','u*_xy','ol*_xy',
            'pt1*_xy','pt_io*_xy','thermal_driving',
            'sa1*_xy','sa_io*_xy','haline_driving'
            ]#,'usws*_xy','vsws*_xy']
plotvar_zlevel = ['pt','sa','u','v']
#tplot = np.arange(40,99,5)
#tplot = np.arange(5,15,5)

#palm.plot_tseries([filedir1], [runname], plotvar_t)
if plot_cross:
    #palm.plot_tseries_cross(diri, run, runlabel, tav=tav_ts, teval=[tplot,tplot],
    #                        plotvar = ['geoice','Umax_i'], col = colorVal, overwrite=True)
    palm.plot_tseries_cross(diri, run, runlabel, tav=tav_ts, teval=[tplot,tplot],
                            plotvar = ['thermal_driving','gamma_T_2m'], col = colorVal, overwrite=True)
if plot_profiles:
    #palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [tplot[0]], ops=norm_pr,zlim=[zmax,0])
    palm.plot_pr(diri, run, plotvar_flux,
                 runlabel=runlabel,  
                 teval = [tprofile,tprofile], 
                 ops=norm_pr, tav = tav_pr, 
                 #xscale_input = xscale_input, xscale_label = xscale_label, zscale = 'Ekman', 
                 zlim=[zmax,0], col=colorVal)#, xlim = [0,20])
    #palm.plot_pr([filedir1],[runname],plotvar = plotvar_pr,teval = [min(tplot),max(tplot)],tall=True, ops=norm_pr,zlim=[zmax,0])
#palm.plot_tseries_zlevel([filedir1], [runname], plotvar_zlevel,norm=['' for i in plotvar_zlevel], zeval=[-10.], tlim=[tplot[0],tplot[-1]])

#palm.plot_pr_varcmp([filedir1],[runname],plotvar_varcmp,teval = tplot,zlim=[zmax,0])

#palm.plot_uv_vector(filedir1,run1,zval = [-2*512/3,0], teval = tplot)
#palm.plot_pr_TS(filedir1,run1,zval = [-2*512/3,0], teval = tplot)
if plot_slices:
    for i,filedir in enumerate(diri):
         palm.plot_3d_slice(filedir,run[i],plotvar,teval = [40],slice_dim = 'y',zval=[zmax,0],yval=[64,64])
#        for z in [-1]:#,-5,-10,-20]:
#        #for z in [-1,-5,-10,-20]:
#            palm.plot_3d_slice(filedir,run[i],plotvar,teval = np.arange(40,49,1),slice_dim = 'z',zval=[z,z])
#palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot,slice_dim = 'y',ops=['' for i in plotvar],yval=[125.,125.],zval=[zmax,0])
#palm.plot_3d_slice(filedir1,runname,plotvar = ['u"'],teval = tplot,ops=[''],slice_dim='y')
#for j in [-1.]:
#    palm.plot_3d_slice(filedir1,runname,plotvar,teval = tplot, slice_dim = 'z', zeval=[j,j],norm=norm_3d)
#palm.plot_2d_xy(filedir1,runname,['melt*_xy'],teval = [10.])
if plot_hovmoller:
    palm.plot_hovmoller(diri,run,['e*'],tlim = [tmin,tmax], runlabel=runlabel,zlim=[zmax,0],overwrite=True, clim=[-7.5,-4.25])
    palm.plot_hovmoller(diri,run,['vel_var_ratio'],tlim = [tmin,tmax], runlabel=runlabel,zlim=[zmax,0],overwrite=True, clim=[0.1,2])
