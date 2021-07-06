#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: cbegeman
"""
import math,os
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as tick
import my_array_mod as ar
import my_string_mod as strop
import convert_NCAR_mod as ncar
import difflib
import gsw
import cmocean as cmo
from plot_param_palm import lw1
#global var,data

# parameters
K0 = 273.15 # temperature (K) at 0 degC
pt_offset = (2/3)*64*(0.01/100)
g = -9.81
xi_N = 0.0052
Ri_crit = 0.2
s_yr = 24.*3600.*365.25
s_hr = 3600.
L_i = 330000. # latent heat of fusion [ J/kg ]
cpw = 4218. # heat capacity of seawater [J/kg deg C ]
omega = 7.29212E-5
f = 2*omega*math.sin(math.pi*70/180)*math.cos(1*math.pi/180)

# list of variables on u grid. TODO check e
uvar = ['xu','yv','zu','zv','zu_xy','zu_3d','zE',
        'Ug_i','Umax_i','Umax','U','U"','u','u"','u*','umax','vmax','v','v"',
        'dEdt','E','E*','dedt','e','e*','gamma_T','gammaT','u*2','v*2',
        'Fshear','Ftrans','Fbuoy','Fbuoy_u','Fbuoy_w','FbuoyT','FbuoyS','diss',
        'w*u*u*:dz','w*p*:dz','N1','N2','S2','Ri','Rf',
        'km','kh','ks','l','e*','w*e*','u2v2','dudz','dvdz','dUdz']

# list of variables on scalar grid
scvar = ['x','y','z','alpha_T','drho_oceandz','drho_gswdz','dprho_gswdz',
         'gamma_T_2m','gamma_T','gammaT','pt(0)', 'pt(z_mo)','thermal_driving_i',
         'pt','pt"','pt*2','sa','sa"','prho','rho_ocean','z_BL','z_SL',
         'rho_gsw','prho_gsw','zu1_xy','sa*2']
# list of variables on w grid
wvar = ['melt','melt*_xy','shf*_xy','shf*_xy_av','sasws*_xy','u*_xy','z0*_xy',
        'dz','zw','zw_xy','zw_3d','k_offset_mcph','Sk','km_eff','kh_eff',
        'w','w"','w*2','zmin_w2','w*3','NORM_ws3','dT',
        'wpt','w*pt*','w"pt"','wsa','w*sa*','wu','wv',
        'wb','w"sa"','w*u*','w"u"','w*v*','w"v"','w"theta"0',
        'z0h*_xy','ol','ol*_xy','pt1*_xy','pt_io*_xy','sa1*_xy','sa_io*_xy',
        'w"U"0','w"u"0','w"v"0','w"pt"0','pt1_t','Sw',
        'haline_driving','usws*_xy','vsws*_xy','freezing',
        'b11','b12','b13','b23','b22','b33','vel_var_ratio']
# list of  variables not defined on a grid
novar = ['time','dt','alpha_surface','thermal_driving_infty','sin_alpha']

# derived variables must be contained within a grid list
dervar = ['alpha_surface','sin_alpha','dedt','dudz','dvdz','dUdz','drho_oceandz','drho_gswdz','dprho_gswdz',
          'U','thermal_driving_i','thermal_driving_infty','haline_driving','freezing','zeta','z_BL','z_SL',
          'w"U"0','wb','K_surf','dT','gamma_T','gammaT','gamma_T_2m',
          'Fshear','Ftrans','FbuoyT','FbuoyS','Fbuoy_u','Fbuoy_w','Fbuoy','diss','dEdt',
          'N1','N2','S2','Ri','Rf','pt1_t','Sw','zmin_w2','u2v2','rho_gsw','prho_gsw',
          'b11','b12','b13','b23','b22','b33','vel_var_ratio','km_eff','kh_eff',
          'U"','Ug_i','Umax_i','Umax','u"','v"','w"','pt"','sa"','u2v2','zE']

varlist = uvar + wvar + scvar + novar
varunits = varlist.copy()
vartype = varlist.copy() 
varname = varlist.copy()
varlabel = varlist.copy()
varcmap = ['cmo.speed' for i in varlist]
varscale = ['linear' for i in varlist]
varscale[varlist.index('E*')] = 'log'
varscale[varlist.index('w*2')] = 'log'
varscale[varlist.index('e*')] = 'log'
varscale[varlist.index('melt')] = 'log'
#varscale[varlist.index('tke_all')] = 'symlog'
#varscale[varlist.index('Fshear')] = 'log'
#varscale[varlist.index('Ftrans')] = 'symlog'
#varscale[varlist.index('Fbuoy')] = 'symlog'
#varscale[varlist.index('Fbuoy_u')] = 'symlog'
#varscale[varlist.index('Fbuoy_w')] = 'symlog'

for i in uvar:
    vartype[varlist.index(i)] = 'u'
for i in wvar:
    vartype[varlist.index(i)] = 'w'
for i in scvar:
    vartype[varlist.index(i)] = 'sc'
    
varname[varlist.index('E')] = 'e'
varname[varlist.index('E*')] = 'e_res'
varname[varlist.index('u*')] = 'us'
varname[varlist.index('e*')] = 'e_res'
varname[varlist.index('melt*_xy')] = 'melt'
varname[varlist.index('ol*_xy')] = 'ol'
varname[varlist.index('pt(0)')] = 'pt0'
varname[varlist.index('pt(z_mo)')] = 'pt_zmo'
varname[varlist.index('pt1_t')] = 'pt_mcphee'
varname[varlist.index('pt1*_xy')] = 'pt_mcphee'
varname[varlist.index('pt_io*_xy')] = 'pt_io'
varname[varlist.index('pt*2')] = 'pt2'
varname[varlist.index('sa1*_xy')] = 'sa_mcphee'
varname[varlist.index('sa_io*_xy')] = 'sa_io'
varname[varlist.index('sa*2')] = 'sa2'
varname[varlist.index('shf*_xy')] = 'shf' 
varname[varlist.index('u*2')] = 'u2' 
varname[varlist.index('u*_xy')] = 'us$'
varname[varlist.index('usws*_xy')] = 'uw' 
varname[varlist.index('v*2')] = 'v2' 
varname[varlist.index('vsws*_xy')] = 'vw' 
varname[varlist.index('w*2')] = 'w2' 
varname[varlist.index('w*p*:dz')] = 'F_Ptrans'
varname[varlist.index('w*pt*')] = 'wpt_res' 
varname[varlist.index('w"pt"')] = 'wpt_sgs' 
varname[varlist.index('w"pt"0')] = 'wpt_sgs_surf' 
varname[varlist.index('w*sa*')] = 'wsa_res' 
varname[varlist.index('w"sa"')] = 'wsa_sgs' 
varname[varlist.index('w*u*u*:dz')] = 'F_Utrans'
varname[varlist.index('w"u"0')] = 'uw_sgs_surf' 
varname[varlist.index('w"u"')] = 'uw_sgs' 
varname[varlist.index('w*u*')] = 'uw_res' 
varname[varlist.index('w"v"0')] = 'vw_sgs_surf' 
varname[varlist.index('w"v"')] = 'vw_sgs'
varname[varlist.index('w*v*')] = 'vw_res'

varlabel[varlist.index('alpha_surface')]= r'\alpha'
varlabel[varlist.index('sin_alpha')]    = r'\sin\alpha'
varlabel[varlist.index('b11')]          = r'b_{11}'
varlabel[varlist.index('b12')]          = r'b_{12}'
varlabel[varlist.index('b13')]          = r'b_{13}'
varlabel[varlist.index('b23')]          = r'b_{23}'
varlabel[varlist.index('b22')]          = r'b_{22}'
varlabel[varlist.index('b33')]          = r'b_{33}'
varlabel[varlist.index('diss')]         = r'\vareps'
varlabel[varlist.index('drho_oceandz')] = r'd\rho/dz'
varlabel[varlist.index('drho_gswdz')]   = r'd\rho/dz'
varlabel[varlist.index('dprho_gswdz')]  = r'd\sigma/dz'
varlabel[varlist.index('dt')]           = r'\delta t'
varlabel[varlist.index('e')]            = r'TKE\:SGS'
#varlabel[varlist.index('e*')]           = r'TKE\:\textrm{resolved}'
varlabel[varlist.index('E*')]           = r'\textrm{TKE}'
varlabel[varlist.index('e*')]           = r'\textrm{TKE}'
varlabel[varlist.index('gamma_T_2m')]   = r'\Gamma_{T,der}'
varlabel[varlist.index('gamma_T')]      = r'\Gamma_{T}'
varlabel[varlist.index('gammaT')]       = r'\gamma_{T}'
varlabel[varlist.index('Fbuoy')]        = r'F_{buoy}'
varlabel[varlist.index('Fbuoy_u')]      = r'F_{buoy,u}'
varlabel[varlist.index('Fbuoy_w')]      = r'F_{buoy,w}'
varlabel[varlist.index('FbuoyT')]       = r'F_{buoy,\theta}'
varlabel[varlist.index('FbuoyS')]       = r'F_{buoy,S}'
varlabel[varlist.index('Fshear')]       = r'F_{shear}'
varlabel[varlist.index('Ftrans')]       = r'F_{trans}'
varlabel[varlist.index('haline_driving')]=r'S - S_i'
varlabel[varlist.index('kh_eff')]       = r'k_h'
varlabel[varlist.index('km_eff')]       = r'k_m'
varlabel[varlist.index('kh')]           = r'k_h SGS'
varlabel[varlist.index('ks')]           = r'k_s SGS'
varlabel[varlist.index('km')]           = r'k_m SGS'
varlabel[varlist.index('k_offset_mcph')]= r'k_{off}'
varlabel[varlist.index('melt')]         = r'\textrm{melt rate}'
varlabel[varlist.index('melt*_xy')]     = r'melt \: rate'
varlabel[varlist.index('N2')]           = r'N^2'
varlabel[varlist.index('N1')]           = r'1/N'
varlabel[varlist.index('ol')]           = r'L_O'
varlabel[varlist.index('ol*_xy')]       = r'L_O'
varlabel[varlist.index('prho')]         = r'\sigma'
varlabel[varlist.index('prho_gsw')]     = r'\sigma'
varlabel[varlist.index('pt(0)')]        = 'pt(0)'
varlabel[varlist.index('pt(z_mo)')]     = 'pt(z_{mo})'
varlabel[varlist.index('pt')]           = r'\theta'
varlabel[varlist.index('pt"')]          = r'\theta^\prime'
varlabel[varlist.index('pt1_t')]        = r'\theta_{m}'
varlabel[varlist.index('pt1*_xy')]      = r'\theta_{m}'
varlabel[varlist.index('pt_io*_xy')]    = r'\theta_{b}'
varlabel[varlist.index('pt*2')]         = r'\overline{{\theta^\prime}^2}/\Theta_0^2'
varlabel[varlist.index('rho_ocean')]    = r'\rho'
varlabel[varlist.index('sa')]           = r'S'
varlabel[varlist.index('sa"')]          = r'S^\prime'
varlabel[varlist.index('sa1*_xy')]      = r'S_{\infty}'
varlabel[varlist.index('sa_io*_xy')]    = r'S_{b}'
varlabel[varlist.index('sasws*_xy')]    = r'\overline{w^\prime S^\prime}'
varlabel[varlist.index('shf*_xy')]      = r'\overline{w^\prime \theta^\prime}'
varlabel[varlist.index('Sw')]           = r'S_w'
varlabel[varlist.index('thermal_driving_i')] = r'\Delta \theta_i'
#varlabel[varlist.index('thermal_driving_infty')] = r'T_{\infty} - T_f(S_{\infty})'
varlabel[varlist.index('thermal_driving_infty')] = r'\Delta \theta_{\infty}'
varlabel[varlist.index('time')]         = r't'
varlabel[varlist.index('Ug_i')]         = r'\overline{u}_g'
varlabel[varlist.index('Umax')]         = r'max(\overline{u})'
varlabel[varlist.index('Umax_i')]       = r'max(\overline{u})'
varlabel[varlist.index('u')]            = r'u'
varlabel[varlist.index('u"')]           = r'u^\prime'
varlabel[varlist.index('U"')]           = r'U^\prime'
varlabel[varlist.index('u*')]           = r'u_*'
varlabel[varlist.index('u*2')]          = r'\overline{u^{\prime 2}}'
varlabel[varlist.index('u*_xy')]        = r'u_*'
varlabel[varlist.index('usws*_xy')]     = r'\overline{u^\prime w^\prime}'
varlabel[varlist.index('v')]            = r'v'
varlabel[varlist.index('v"')]           = r'v^\prime'
varlabel[varlist.index('v*2')]          = r'\overline{v^{\prime 2}}'
varlabel[varlist.index('vsws*_xy')]     = r'\overline{v^\prime w^\prime}'
varlabel[varlist.index('vel_var_ratio')]= r'2\overline{w^{\prime 2}}/(\overline{u^{\prime 2}}+\overline{v^{\prime 2}})'
varlabel[varlist.index('w')]            = r'w'
varlabel[varlist.index('w"')]           = r'w^\prime'
varlabel[varlist.index('w*2')]          = r'\overline{w^{\prime 2}}'
varlabel[varlist.index('w*p*:dz')]      = r'F_{Ptrans}'#Transport of resolved TKE due to pressure fluctuations (term in resolved TKE budget)
varlabel[varlist.index('wpt')]          = r'\overline{w^\prime \theta^\prime}'
varlabel[varlist.index('w*pt*')]        = r'\overline{w^\prime \theta^\prime} resolved'
varlabel[varlist.index('w"pt"')]        = r'\overline{w^\prime \theta^\prime} SGS'
varlabel[varlist.index('w*sa*')]        = r'\overline{w^\prime S^\prime} resolved'
varlabel[varlist.index('w"sa"')]        = r'\overline{w^\prime S^\prime} SGS'
varlabel[varlist.index('wsa')]          = r'\overline{w^\prime S^\prime}'
varlabel[varlist.index('wu')]           = r'\overline{u^\prime w^\prime}'
varlabel[varlist.index('w*u*u*:dz')]    = r'F_{Utrans}'#Transport of resolved TKE due to turbulance (term in resolved TKE budget)
varlabel[varlist.index('w"u"')]         = r'\overline{u^\prime w^\prime} SGS'
varlabel[varlist.index('w*u*')]         = r'\overline{u^\prime w^\prime} resolved'
varlabel[varlist.index('wv')]           = r'\overline{v^\prime w^\prime}'
varlabel[varlist.index('w"v"')]         = r'\overline{v^\prime w^\prime} SGS'
varlabel[varlist.index('w*v*')]         = r'\overline{v^\prime w^\prime} resolved'
varlabel[varlist.index('x')]            = r'x'
varlabel[varlist.index('xu')]           = r'x'
varlabel[varlist.index('y')]            = r'y'
varlabel[varlist.index('yv')]           = r'y'
varlabel[varlist.index('z_BL')]         = r'z_{BL}'
varlabel[varlist.index('z')]            = r'z'
varlabel[varlist.index('zu')]           = r'z'
varlabel[varlist.index('zu_xy')]        = r'z'
varlabel[varlist.index('zu_3d')]        = r'z'
varlabel[varlist.index('zw')]           = r'z'

varunits[varlist.index('alpha_surface')]= r'^{\circ}'
varunits[varlist.index('sin_alpha')]    = r''
varunits[varlist.index('b11')]          = r''
varunits[varlist.index('b12')]          = r''
varunits[varlist.index('b13')]          = r''
varunits[varlist.index('b23')]          = r''
varunits[varlist.index('b22')]          = r''
varunits[varlist.index('gammaT')]       = r''
varunits[varlist.index('gamma_T')]      = r''
varunits[varlist.index('gamma_T_2m')]   = r''
varunits[varlist.index('b33')]          = r''
varunits[varlist.index('diss')]         = r'm^2/s^3'
varunits[varlist.index('dz')]           = r'm^{-1}'
varunits[varlist.index('dt')]           = r's'
varunits[varlist.index('e')]            = r'm^2\:s^{-2}'
varunits[varlist.index('E*')]           = r'm^2\:s^{-2}'
varunits[varlist.index('e*')]           = r'm^2\:s^{-2}'
varunits[varlist.index('Fbuoy')]        = r'm^2\:s^{-3}'
varunits[varlist.index('Fbuoy_u')]      = r'm^2\:s^{-3}'
varunits[varlist.index('Fbuoy_w')]      = r'm^2\:s^{-3}'
varunits[varlist.index('FbuoyT')]       = r'm^2\:s^{-3}'
varunits[varlist.index('FbuoyS')]       = r'm^2\:s^{-3}'
varunits[varlist.index('Fshear')]       = r'm^2\:s^{-3}'
varunits[varlist.index('Ftrans')]       = r'm^2\:s^{-3}'
varunits[varlist.index('haline_driving')]=r'psu'
varunits[varlist.index('kh_eff')]       = r'm^2\:s^{-1}'
varunits[varlist.index('km_eff')]       = r'm^2\:s^{-1}'
varunits[varlist.index('kh')]           = r'm^2\:s^{-1}'
varunits[varlist.index('ks')]           = r'm^2\:s^{-1}'
varunits[varlist.index('km')]           = r'm^2\:s^{-1}'
varunits[varlist.index('k_offset_mcph')]= r''
varunits[varlist.index('l')]            = r'm'
varunits[varlist.index('melt')]         = r'm/s'
varunits[varlist.index('melt*_xy')]     = r'm/s'
varunits[varlist.index('N2')]           = r's^{-2}'
varunits[varlist.index('N1')]           = r'min'
varunits[varlist.index('ol')]           = r'm'
varunits[varlist.index('ol*_xy')]       = r'm'
varunits[varlist.index('prho')]         = r'kg m^{-3}'
varunits[varlist.index('prho_gsw')]     = r'kg m^{-3}'
varunits[varlist.index('pt')]           = r'^{\circ}C'
varunits[varlist.index('pt"')]          = r'^{\circ}C'
varunits[varlist.index('pt1_t')]        = r'^{\circ}C'
varunits[varlist.index('pt1*_xy')]      = r'^{\circ}C'
varunits[varlist.index('pt_io*_xy')]    = r'^{\circ}C'
varunits[varlist.index('pt*2')]         = r'^{\circ}C^2'
varunits[varlist.index('rho_ocean')]    = r'kg m^{-3}'
varunits[varlist.index('sa')]           = r'psu'
varunits[varlist.index('sa"')]          = r'psu'
varunits[varlist.index('sa1*_xy')]      = r'psu'
varunits[varlist.index('sa_io*_xy')]    = r'psu'
varunits[varlist.index('sasws*_xy')]    = r'psu m s^{-1}'
varunits[varlist.index('shf*_xy')]      = r'K m s^{-1}'
varunits[varlist.index('Sw')]           = r'(m^3 s^{-2})^{-0.5}'
varunits[varlist.index('thermal_driving_i')]=r'^{\circ}C'
varunits[varlist.index('thermal_driving_infty')]=r'^{\circ}C'
varunits[varlist.index('time')]         = r's'
varunits[varlist.index('u')]            = r'm s^{-1}'
varunits[varlist.index('u"')]           = r'm s^{-1}'
varunits[varlist.index('U"')]           = r'm s^{-1}'
varunits[varlist.index('U')]            = r'm s^{-1}'
varunits[varlist.index('Ug_i')]         = r'm s^{-1}'
varunits[varlist.index('Umax')]         = r'm s^{-1}'
varunits[varlist.index('Umax_i')]         = r'm s^{-1}'
varunits[varlist.index('u*')]           = r'm\:s^{-1}'
varunits[varlist.index('u*_xy')]        = r'm s^{-1}'
varunits[varlist.index('u*2')]          = r'm^2 s^{-2}'
varunits[varlist.index('umax')]         = r'm s^{-1}'
varunits[varlist.index('usws*_xy')]     = r'm^2 s^{-2}'
varunits[varlist.index('v')]            = r'm s^{-1}'
varunits[varlist.index('v"')]           = r'm s^{-1}'
varunits[varlist.index('vel_var_ratio')]= r''
varunits[varlist.index('v*2')]          = r'm^2 s^{-2}'
varunits[varlist.index('vmax')]         = r'm s^{-1}'
varunits[varlist.index('vsws*_xy')]     = r'm^2 s^{-2}'
varunits[varlist.index('w')]            = r'm\:s^{-1}'
varunits[varlist.index('w"')]           = r'm\:s^{-1}'
varunits[varlist.index('w*2')]          = r'm^2 s^{-2}'
varunits[varlist.index('w*3')]          = r'm^2 s^{-2}'
varunits[varlist.index('w*p*:dz')]      = r'm^2 s^{-2}'
varunits[varlist.index('wpt')]          = r'K m s^{-1}'
varunits[varlist.index('w"pt"')]        = r'K m s^{-1}'
varunits[varlist.index('w*pt*')]        = r'K m s^{-1}'
varunits[varlist.index('w"pt"0')]        = r'K m s^{-1}'
varunits[varlist.index('w*sa*')]        = r'psu\:m\:s^{-1}'
varunits[varlist.index('w"sa"')]        = r'psu\:m\:s^{-1}'
varunits[varlist.index('wsa')]          = r'psu\:m\:s^{-1}'
varunits[varlist.index('w*u*u*:dz')]    = r'm^2\:s^{-2}'
varunits[varlist.index('w"u"')]         = r'm^2\:s^{-2}'
varunits[varlist.index('wu')]           = r'm^2\:s^{-2}'
varunits[varlist.index('wv')]           = r'm^2\:s^{-2}'
varunits[varlist.index('w"v"')]         = r'm^2\:s^{-2}'
varunits[varlist.index('x')]            = r'm'
varunits[varlist.index('xu')]           = r'm'
varunits[varlist.index('y')]            = r'm'
varunits[varlist.index('yv')]           = r'm'
varunits[varlist.index('z_BL')]         = r'm'
varunits[varlist.index('z')]            = r'm'
varunits[varlist.index('zu')]           = r'm'
varunits[varlist.index('zu_xy')]        = r'm'
varunits[varlist.index('zu_3d')]        = r'm'
varunits[varlist.index('zv')]           = r'm'
varunits[varlist.index('zw')]           = r'm'
varunits[varlist.index('zw_xy')]        = r'm'
varunits[varlist.index('zw_3d')]        = r'm'
varunits[varlist.index('z_SL')]         = r'm'
#units derived from other terms here
varunits[varlist.index('drho_oceandz')] = varunits[varlist.index('rho_ocean')]+varunits[varlist.index('dz')] 
varunits[varlist.index('drho_gswdz')]   = varunits[varlist.index('rho_ocean')]+varunits[varlist.index('dz')] 
varunits[varlist.index('dprho_gswdz')]  = varunits[varlist.index('prho')]+varunits[varlist.index('dz')] 

#varcmap[:] = 'cmo.speed'
varcmap[varlist.index('shf*_xy')] = 'cmo.thermal'
varcmap[varlist.index('sasws*_xy')] = 'cmo.haline'
#varcmap[varlist.index('thermal_driving')] = 'cmo.thermal'
varcmap[varlist.index('haline_driving')] = 'cmo.haline'
varcmap[varlist.index('pt')] = 'cmo.thermal'
varcmap[varlist.index('pt"')] = 'cmo.thermal'
varcmap[varlist.index('pt1*_xy')] = 'cmo.thermal'
varcmap[varlist.index('pt_io*_xy')] = 'cmo.thermal'
varcmap[varlist.index('sa')] = 'cmo.haline'
varcmap[varlist.index('sa"')] = 'cmo.haline'
varcmap[varlist.index('sa1*_xy')] = 'cmo.haline'
varcmap[varlist.index('sa_io*_xy')] = 'cmo.haline'
varcmap[varlist.index('rho_ocean')] = 'cmo.dense'
varcmap[varlist.index('prho')] = 'cmo.dense'
varcmap[varlist.index('u*_xy')] = 'cmo.speed'
varcmap[varlist.index('usws*_xy')] = 'cmo.balance'
varcmap[varlist.index('vsws*_xy')] = 'cmo.balance'
varcmap[varlist.index('U')] = 'cmo.speed'
varcmap[varlist.index('u')] = 'cmo.balance'
varcmap[varlist.index('v')] = 'cmo.balance'
varcmap[varlist.index('w')] = 'cmo.balance'
varcmap[varlist.index('u"')] = 'cmo.balance'
varcmap[varlist.index('v"')] = 'cmo.balance'
varcmap[varlist.index('w"')] = 'cmo.balance'
varcmap[varlist.index('melt*_xy')] = 'cmo.balance'
varcmap[varlist.index('ol*_xy')] = 'cmo.speed'

varsname = (['variance','velocity',
             #'vel_var_ratio','gamma_T','gammaT','Ug_i','Umax','Umax_i',
             #'dT','thermal_driving','haline_driving',
             'k_all','hor_vert_variance',
             #'w"U"0','S2','N2','u2v2',
             'momflux_u','momflux_v','momflux_z',
             'heatflux_z','saltflux_z',
             #'dEdt','dedt','dudz','dvdz','dUdz','drho_oceandz',
             #'drho_gswdz','dprho_gswdz',
             'Fbuoy_uw','tke','tke_all'])
             #'Fshear','Fbuoy','FbuoyT','FbuoyS'])
varsvars = [list() for i in varsname]
vars_axis_label = ['' for i in varsname]
vars_ls = [list() for i in varsname]
vars_lw = [list() for i in varsname]
varsvars[varsname.index('variance')] = ['u*2','v*2','w*2']
vars_axis_label[varsname.index('variance')] = r'Velocity \: variance \: (m^2/s^2)'
vars_ls[varsname.index('variance')] = ['-','--',':']
vars_lw[varsname.index('variance')] = [lw1 for i in varsvars[varsname.index('variance')]]
varsvars[varsname.index('hor_vert_variance')] = ['u2v2','w*2']
vars_axis_label[varsname.index('hor_vert_variance')] = r'Velocity \: variance \: (m^2/s^2)'
vars_ls[varsname.index('hor_vert_variance')] = ['-','--']
vars_lw[varsname.index('hor_vert_variance')] = [lw1 for i in varsvars[varsname.index('hor_vert_variance')]]
varsvars[varsname.index('momflux_u')] = ['wu','w*u*','w"u"']
vars_axis_label[varsname.index('momflux_u')] = r'uw \: (m^2/s^2)'
vars_ls[varsname.index('momflux_u')] = ['-','--',':']
vars_lw[varsname.index('momflux_u')] = [lw1 for i in varsvars[varsname.index('momflux_u')]]
varsvars[varsname.index('momflux_v')] = ['wv','w*v*','w"v"']
vars_axis_label[varsname.index('momflux_v')] = r'vw \: (m^2/s^2)'
vars_ls[varsname.index('momflux_v')] = ['-','--',':']
vars_lw[varsname.index('momflux_v')] = [lw1 for i in varsvars[varsname.index('momflux_v')]]
varsvars[varsname.index('momflux_z')] = ['wv','wu']
vars_axis_label[varsname.index('momflux_z')] = r'\overline{u_i^\prime w^\prime} \: (m^2\:s^{-2})'
vars_ls[varsname.index('momflux_z')] = ['-','--']
vars_lw[varsname.index('momflux_z')] = [lw1,lw1]
varsvars[varsname.index('heatflux_z')] = ['wpt','w*pt*','w"pt"']
vars_axis_label[varsname.index('heatflux_z')] = varlabel[varname.index('wpt')] + '\: (' + varunits[varname.index('wpt')] + ')' 
vars_ls[varsname.index('heatflux_z')] = ['-','--',':']
vars_lw[varsname.index('heatflux_z')] = [lw1 for i in varsvars[varsname.index('heatflux_z')]]
varsvars[varsname.index('saltflux_z')] = ['wsa','w*sa*','w"sa"']
vars_axis_label[varsname.index('saltflux_z')] = varlabel[varname.index('wsa')] + '\: (' + varunits[varname.index('wsa')] + ')' 
vars_ls[varsname.index('saltflux_z')] = ['-','--',':']
vars_lw[varsname.index('saltflux_z')] = [lw1 for i in varsvars[varsname.index('saltflux_z')]]
varsvars[varsname.index('tke')] = ['Fshear','Fbuoy']
vars_axis_label[varsname.index('tke')] = r'dE/dt \: (m^2/s^3)' 
vars_ls[varsname.index('tke')] = ['-','--']
vars_lw[varsname.index('tke')] = [lw1 for i in varsvars[varsname.index('tke')]]

varsvars[varsname.index(       'k_all')] = ['km','kh','ks']
vars_axis_label[varsname.index('k_all')] = r'k \: (m^2/s)'
vars_ls[varsname.index(        'k_all')] = ['-','--',':']
vars_lw[varsname.index(        'k_all')] = [lw1,lw1,lw1]

#varsvars[varsname.index('tke_all')] = ['Fshear','Fbuoy','FbuoyT','FbuoyS','w*p*:dz','w*u*u*:dz']
varsvars[varsname.index('tke_all')] = ['Fshear','Fbuoy','Ftrans']
vars_axis_label[varsname.index('tke_all')] = r'dE/dt \: (m^2/s^3)' 
vars_ls[varsname.index('tke_all')] = ['-','--',':']
vars_lw[varsname.index('tke_all')] = [lw1,lw1,lw1]

varsvars[varsname.index('Fbuoy_uw')] = ['Fbuoy','Fbuoy_w','Fbuoy_u']
vars_axis_label[varsname.index('Fbuoy_uw')] = r'F_{buoy} \: (m^2\:s^{-3})' 
vars_ls[varsname.index('Fbuoy_uw')] = ['-','--',':']
vars_lw[varsname.index('Fbuoy_uw')] = [lw1,lw1,lw1]

varsvars[varsname.index('velocity')] = ['u','v']
vars_axis_label[varsname.index('velocity')] = r'u_i \: (m/s)'
vars_ls[varsname.index('velocity')] = ['-','--']
vars_lw[varsname.index('velocity')] = [lw1,lw1]

# include derived variable proxies so that shape is extracted in slice_var
source_var = varlist.copy()
source_var[varlist.index('w"U"0')]           = 'w"u"0'
source_var[varlist.index('vel_var_ratio')]   = 'u*2'
source_var[varlist.index('u2v2')]            = 'u*2'
source_var[varlist.index('thermal_driving_i')] = 'u*'
source_var[varlist.index('thermal_driving_infty')] = 'u*'
source_var[varlist.index('dT')]              = 'melt'
source_var[varlist.index('haline_driving')]  = 'sa_io*_xy'
source_var[varlist.index('Fshear')]          = 'u'
source_var[varlist.index('Ftrans')]          = 'u'
source_var[varlist.index('Fbuoy')]           = 'wpt'
source_var[varlist.index('Fbuoy_u')]         = 'u"pt"'
source_var[varlist.index('Fbuoy_w')]         = 'wpt'
source_var[varlist.index('FbuoyT')]          = 'wpt'
source_var[varlist.index('FbuoyS')]          = 'wpt'
source_var[varlist.index('Ug_i')]            = 'u'
source_var[varlist.index('Umax_i')]          = 'u'
source_var[varlist.index('Umax')]            = 'u'
source_var[varlist.index('S2')]              = 'u'
source_var[varlist.index('Ri')]              = 'u'
source_var[varlist.index('Rf')]              = 'u'
source_var[varlist.index('N2')]              = 'rho_ocean'
source_var[varlist.index('gamma_T_2m')]      = 'u'
source_var[varlist.index('gammaT')]          = 'u*'
source_var[varlist.index('gamma_T')]         = 'u*'
source_var[varlist.index('dEdt')]            = 'e'
source_var[varlist.index('dedt')]            = 'e'
source_var[varlist.index('dudz')]            = 'u'
source_var[varlist.index('dvdz')]            = 'u'
source_var[varlist.index('dUdz')]            = 'u'
source_var[varlist.index('U')]              = 'u'
source_var[varlist.index('U"')]              = 'u'
source_var[varlist.index('drho_oceandz')]    = 'rho_ocean'
source_var[varlist.index('drho_gswdz')]      = 'pt'
source_var[varlist.index('dprho_gswdz')]     = 'pt'
source_var[varlist.index('km_eff')]     = 'km'


datafilelist = ['parameter','x','y','z','ts','pr','2d','3d'] # corresponds to data_type
filename = ['','','','',
            'DATA_1D_TS_NETCDF',
            'DATA_1D_PR_NETCDF',
            'DATA_2D_XY_NETCDF',
            'DATA_3D_NETCDF']# add averaging options
data_dim = [[''],['x'],['y'],['z'],['t'],['t','z'],['t','z','x','y'],['t','z','x','y']]


