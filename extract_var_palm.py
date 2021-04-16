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

#global var,data

# parameters
K0 = 273.15 # temperature (K) at 0 degC
g = -9.81
xi_N = 0.0052
Ri_crit = 0.2
s_yr = 24.*3600.*365.25
s_hr = 3600.

# plotting parameters
lw1 = 2 
lw2 = 3
lw3 = 1 

# list of variables on u grid. TODO check e
uvar = ['xu','yv','zu','zv','zu_xy','zu_3d',
        'U','u','u"','u*','umax','vmax','v','v"',
        'dEdt','E','E*','dedt','e','e*','u*2','v*2',
        'Fshear','Ftrans','Fbuoy','FbuoyT','FbuoyS','diss',
        'w*u*u*:dz','w*p*:dz','N1','N2','S2','Ri','Rf',
        'km','kh','l','e*','w*e*','u2v2','dudz','dvdz','dUdz']

# list of variables on scalar grid
scvar = ['x','y','z','alpha_T','drho_oceandz','drho_gswdz','dprho_gswdz',
         'pt','pt"','pt*2','sa','sa"','prho','rho_ocean','z_BL','z_SL',
         'rho_gsw','prho_gsw','zu1_xy','sa*2']
# list of variables on w grid
wvar = ['melt','melt*_xy','shf*_xy','shf*_xy_av','sasws*_xy','u*_xy','z0*_xy',
        'dz','zw','zw_xy','zw_3d','k_offset_mcph','Sk','km_eff','kh_eff',
        'w','w"','w*2','zmin_w2','w*3','NORM_ws3','dT',
        'wpt','w*pt*','w"pt"','wsa','w*sa*','wu','wv',
        'wb','w"sa"','w*u*','w"u"','w*v*','w"v"','w"theta"0',
        'z0h*_xy','ol','ol*_xy','pt1*_xy','pt_io*_xy','sa1*_xy','sa_io*_xy',
        'w"U"0','w"u"0','w"v"0','w"pt"0','thermal_driving','pt1_t','Sw',
        'haline_driving','usws*_xy','vsws*_xy','freezing',
        'b11','b12','b13','b23','b22','b33','vel_var_ratio']
# list of  variables not defined on a grid
novar = ['time','dt']

# derived variables must be contained within a grid list
dervar = ['dedt','dudz','dvdz','dUdz','drho_oceandz','drho_gswdz','dprho_gswdz',
          'U','thermal_driving','haline_driving','freezing','zeta','z_BL','z_SL',
          'w"U"0','wb','K_surf','dT',
          'Fshear','Ftrans','FbuoyT','FbuoyS','Fbuoy','diss','dEdt',
          'N1','N2','S2','Ri','Rf','pt1_t','Sw','zmin_w2','u2v2','rho_gsw','prho_gsw',
          'b11','b12','b13','b23','b22','b33','vel_var_ratio','km_eff','kh_eff',
          'u"','v"','w"','pt"','sa"','u2v2']

varlist = uvar + wvar + scvar + novar
varunits = varlist.copy()
vartype = varlist.copy() 
varname = varlist.copy()
varlabel = varlist.copy()
varcmap = ['cmo.speed' for i in varlist]

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
varlabel[varlist.index('e')]            = r'TKE SGS'
varlabel[varlist.index('e*')]           = r'TKE resolved'
varlabel[varlist.index('Fbuoy')]        = r'F_{buoy}'
varlabel[varlist.index('FbuoyT')]       = r'F_{buoy,\theta}'
varlabel[varlist.index('FbuoyS')]       = r'F_{buoy,S}'
varlabel[varlist.index('Fshear')]       = r'F_{shear}'
varlabel[varlist.index('Ftrans')]       = r'F_{trans}'
varlabel[varlist.index('haline_driving')]=r'S - S_i'
varlabel[varlist.index('kh_eff')]       = r'k_h'
varlabel[varlist.index('km_eff')]       = r'k_m'
varlabel[varlist.index('kh')]           = r'k_h SGS'
varlabel[varlist.index('km')]           = r'k_m SGS'
varlabel[varlist.index('k_offset_mcph')]= r'k_{off}'
varlabel[varlist.index('melt')]         = r'melt \: rate'
varlabel[varlist.index('melt*_xy')]     = r'melt \: rate'
varlabel[varlist.index('N2')]           = r'N^2'
varlabel[varlist.index('N1')]           = r'1/N'
varlabel[varlist.index('ol')]           = r'L_O'
varlabel[varlist.index('ol*_xy')]       = r'L_O'
varlabel[varlist.index('prho')]         = r'\sigma'
varlabel[varlist.index('prho_gsw')]     = r'\sigma'
varlabel[varlist.index('pt')]           = r'\theta'
varlabel[varlist.index('pt"')]          = r'\theta^\prime'
varlabel[varlist.index('pt1_t')]        = r'\theta_{m}'
varlabel[varlist.index('pt1*_xy')]      = r'\theta_{m}'
varlabel[varlist.index('pt_io*_xy')]    = r'\theta_{b}'
varlabel[varlist.index('pt*2')]         = r'\overline{{\theta^\prime}^2}/\Theta_0^2'
varlabel[varlist.index('rho_ocean')]    = r'\rho'
varlabel[varlist.index('sa')]           = r'S'
varlabel[varlist.index('sa"')]          = r'S^\prime'
varlabel[varlist.index('sa1*_xy')]      = r'S_{infty}'
varlabel[varlist.index('sa_io*_xy')]    = r'S_{b}'
varlabel[varlist.index('sasws*_xy')]    = r'\overline{w^\prime S^\prime}'
varlabel[varlist.index('shf*_xy')]      = r'\overline{w^\prime \theta^\prime}'
varlabel[varlist.index('Sw')]           = r'S_w'
varlabel[varlist.index('thermal_driving')] = r'T - T_f(S_i)'
varlabel[varlist.index('time')]         = r't'
varlabel[varlist.index('u')]            = r'u'
varlabel[varlist.index('u"')]           = r'u^\prime'
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
varlabel[varlist.index('z_BL')]         = r'z_{BL}$'
varlabel[varlist.index('z')]            = r'z'
varlabel[varlist.index('zu')]           = r'z'
varlabel[varlist.index('zu_xy')]        = r'z'
varlabel[varlist.index('zu_3d')]        = r'z'
varlabel[varlist.index('zw')]           = r'z'

varunits[varlist.index('b11')]          = r''
varunits[varlist.index('b12')]          = r''
varunits[varlist.index('b13')]          = r''
varunits[varlist.index('b23')]          = r''
varunits[varlist.index('b22')]          = r''
varunits[varlist.index('b33')]          = r''
varunits[varlist.index('diss')]         = r'm^2/s^3'
varunits[varlist.index('dz')]           = r'm^{-1}'
varunits[varlist.index('dt')]           = r's'
varunits[varlist.index('e')]            = r'm^2/s^2'
varunits[varlist.index('e*')]           = r'm^2/s^2'
varunits[varlist.index('Fbuoy')]        = r'm^2/s^3'
varunits[varlist.index('FbuoyT')]       = r'm^2/s^3'
varunits[varlist.index('FbuoyS')]       = r'm^2/s^3'
varunits[varlist.index('Fshear')]       = r'm^2/s^3'
varunits[varlist.index('Fshear')]       = r'm^2/s^3'
varunits[varlist.index('haline_driving')]=r'psu'
varunits[varlist.index('kh_eff')]       = r'm^2/s'
varunits[varlist.index('km_eff')]       = r'm^2/s'
varunits[varlist.index('kh')]           = r'm^2/s'
varunits[varlist.index('km')]           = r'm^2/s'
varunits[varlist.index('k_offset_mcph')]= r''
varunits[varlist.index('l')]            = r'm'
varunits[varlist.index('melt')]         = r'm/s'
varunits[varlist.index('melt*_xy')]     = r'm/s'
varunits[varlist.index('N2')]           = r's^{-1}'
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
varunits[varlist.index('thermal_driving')]=r'^{\circ}C'
varunits[varlist.index('time')]         = r's'
varunits[varlist.index('u')]            = r'm s^{-1}'
varunits[varlist.index('u"')]           = r'm s^{-1}'
varunits[varlist.index('U')]            = r'm s^{-1}'
varunits[varlist.index('u*')]           = r'm s^{-1}'
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
varunits[varlist.index('w')]            = r'm s^{-1}/s'
varunits[varlist.index('w"')]           = r'm s^{-1}/s'
varunits[varlist.index('w*2')]          = r'm^2 s^{-2}'
varunits[varlist.index('w*3')]          = r'm^2 s^{-2}'
varunits[varlist.index('w*p*:dz')]      = r'm^2 s^{-2}'
varunits[varlist.index('wpt')]          = r'K m s^{-1}'
varunits[varlist.index('w"pt"')]        = r'K m s^{-1}'
varunits[varlist.index('w*pt*')]        = r'K m s^{-1}'
varunits[varlist.index('w"pt"0')]        = r'K m s^{-1}'
varunits[varlist.index('w*sa*')]        = r'psu m s^{-1}'
varunits[varlist.index('w"sa"')]        = r'psu m s^{-1}'
varunits[varlist.index('wsa')]          = r'psu m s^{-1}'
varunits[varlist.index('w*u*u*:dz')]    = r'm^2 s^{-2}'
varunits[varlist.index('w"u"')]         = r'm^2 s^{-2}'
varunits[varlist.index('wu')]           = r'm^2 s^{-2}'
varunits[varlist.index('wv')]           = r'm^2 s^{-2}'
varunits[varlist.index('w"v"')]         = r'm^2 s^{-2}'
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
varcmap[varlist.index('thermal_driving')] = 'cmo.thermal'
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

varsname = (['variance','vel_var_ratio',
             'dT','thermal_driving','haline_driving',
             'hor_vert_variance','S2','N2','u2v2',
             'w"U"0','momflux_u','momflux_v','momflux_z',
             'heatflux_z','saltflux_z',
             'dEdt','dedt','dudz','dvdz','dUdz','drho_oceandz',
             'drho_gswdz','dprho_gswdz',
             'tke','tke_all','Fshear','Fbuoy','FbuoyT','FbuoyS','velocity'])
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
vars_axis_label[varsname.index('momflux_z')] = r'u_i w \: (m^2/s^2)'
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

varsvars[varsname.index('tke_all')] = ['Fshear','Fbuoy','FbuoyT','FbuoyS','w*p*:dz','w*u*u*:dz']
vars_axis_label[varsname.index('tke_all')] = r'dE/dt \: (m^2/s^3)' 
vars_ls[varsname.index('tke_all')] = ['-','--','--',':','--',':']
vars_lw[varsname.index('tke_all')] = [lw2,lw2,lw1,lw1,lw3,lw3]

varsvars[varsname.index('velocity')] = ['u','v']
vars_axis_label[varsname.index('velocity')] = r'u_i \: (m/s)'
vars_ls[varsname.index('velocity')] = ['-','--']
vars_lw[varsname.index('velocity')] = [lw1,lw1]

# include derived variable proxies so that shape is extracted in slice_var
varsvars[varsname.index('w"U"0')] = ['w"u"0']
varsvars[varsname.index('vel_var_ratio')] = ['u*2']
varsvars[varsname.index('u2v2')] = ['u*2']
varsvars[varsname.index('thermal_driving')] = ['pt_io*_xy']
varsvars[varsname.index('dT')] = ['melt']
varsvars[varsname.index('haline_driving')] = ['sa_io*_xy']
varsvars[varsname.index('Fshear')] = ['u']
varsvars[varsname.index('Fbuoy')] = ['wpt']
varsvars[varsname.index('FbuoyT')] = ['wpt']
varsvars[varsname.index('FbuoyS')] = ['wpt']
varsvars[varsname.index('S2')] = ['u']
varsvars[varsname.index('N2')] = ['rho_ocean']
varsvars[varsname.index('dEdt')] = ['e']
varsvars[varsname.index('dedt')] = ['e']
varsvars[varsname.index('dudz')] = ['u']
varsvars[varsname.index('dvdz')] = ['u']
varsvars[varsname.index('dUdz')] = ['u']
varsvars[varsname.index('drho_oceandz')] = ['rho_ocean']
varsvars[varsname.index('drho_gswdz')] = ['pt']
varsvars[varsname.index('dprho_gswdz')] = ['pt']

datafilelist = ['x','y','z','ts','pr','2d','3d'] # corresponds to data_type
filename = ['','','',
            'DATA_1D_TS_NETCDF',
            'DATA_1D_PR_NETCDF',
            'DATA_2D_XY_NETCDF',
            'DATA_3D_NETCDF']# add averaging options
data_dim = [['x'],['y'],['z'],['t'],['t','z'],['t','z','x','y'],['t','z','x','y']]

def load_data(filedir,data_type='pr',av=False,coupled=False):
    file1 = filename[datafilelist.index(data_type)]
    if av:
        file1 = file1[:-7] + '_AV' + file1[-7:]
    if coupled:
        file1 = file1 + '_O'
    data = Dataset(filedir+file1,'r')
    return data

def grid_var(var,data_type='pr',grid='sc'):
    if (var == 'z'):
        if (grid == 'sc') or (grid == 'u'):
            if data_type == '3d':
                var = 'zu_3d'                
            elif data_type == '2d':
                var = 'zu_xy'
            else:
                var = 'zu'
        else:
            if data_type == '3d':
                var = 'zw_3d'                
            elif data_type == '2d':
                var = 'zw_xy'
            else:
                var = 'zw'
        data_type = 'z'
    elif (var == 'x'):
        if (grid == 'sc') or (grid == 'u'):
            var = 'xu'
        else:
            var= 'x'
        data_type = 'x'
    elif (var== 'y'):
        if (grid == 'sc') or (grid == 'u'):
            var= 'yv'
        else:
            var= 'y'
        data_type = 'y'
    return var,data_type

def extract_var(data,varname,data_type='pr',ops=[],
                grid='',keep='',
                tval=[9999.,9999.],xval=[9999.,9999.],
                yval=[9999.,9999.],zval=[9999.,9999.],
                p0 = 10., filedir=''):

    print('extract_var: '+varname+', data_type= '+data_type+', grid= ',grid)
    if (varname == 'z') or (varname == 'x') or (varname == 'y'):
        grid = vartype[varlist.index(varname)] 
        varname,data_type = grid_var(varname,data_type=data_type,grid=grid)
    
    slice_obj,varaxes = slice_var(data,varname,data_type=data_type,
                        tval=tval,xval=xval,yval=yval,zval=zval,grid=grid)
    if varname in dervar:
        var1 = derived_var(data,varname,slice_obj,filedir=filedir,p0=p0,
                           data_type=data_type)
    else:
        if not(varname in data.variables.keys()):
            print('variable ' + varname + ' not in file')
            return [],''
        var1 = data.variables[varname][slice_obj]
    
    shp = np.shape(var1);
    sq = ()
    for idx,i in enumerate(varaxes):
        if (i not in keep) and (shp[idx] == 1):
            sq += (int(idx),)
    if len(sq) > 0:
        var1 = np.squeeze(var1,axis=sq);
        del sq
    
    # perform operations on var 
    varlabel1 = varlabel[varlist.index(varname)]
    varlabel2 = r''
    
    if varname[0:2] == 'pt':
        var1 = np.subtract(var1,K0)
        varunits[varlist.index(varname)] = r'^{\circ}C'
    elif 'melt' in varname:
        var1 = np.multiply(var1,s_yr)
        axis_label = varlabel[varlist.index(varname)]+ ' (m/yr)'
    if 'shf' in ops:
        var2,varlabel2 = extract_var(data,'shf',data_type='ts')
    elif 'rho0' in ops:
        var2,varlabel2 = extract_var(data,'rho_ocean',
                         tval=[0,0],zval=[0,0],data_type='pr')
    elif ('mean' in ops) and len(var1) > 1:
        # issue with using axis formulation for dimension 1 arrays
        var2 = np.mean(var1);#,axis = var_axes.index(ops[0]))
        varlabel2 = r'\overline{'+varlabel[varlist.index(varname)]+r'}'
        if ('diff' not in ops) and ('norm' not in ops):
           var1 = var2.copy() # if no further operation, output mean
    elif 'std' in ops:
        var2 = np.std(var1,axis=var_axes.index(ops[0]))
        varlabel2 = r'\sigma_{'+varlabel[varlist.index(varname)]+'}'
        if ('diff' not in ops) and ('norm' not in ops):
           var1 = var2.copy()
    elif 'sasws' in ops:
        var2,varlabel2 = extract_var(data,'sasws',data_type='ts')
    elif 'ws3' in ops:
        var2,varlabel2 = extract_var(data,'NORM_ws3',tval=tval,data_type='pr')
    elif 'hr' in ops:
        var1 = np.divide(var1,3600.)
        varunits[varlist.index(varname)] = 'h'
    if 'ptfr' in ops:
        #sa_far = data.variables['sa'][0,0]
        ##p0 = data.variables['hyp'][1,-1]
        #t_io = gsw.t_freezing(sa_far,p0,0)
        #pt_io = gsw.pt0_from_t(sa_far,t_io,p0)
        #var1 = np.subtract(var1,pt_io)
        #axis_label = varlabel[varlist.index(var)]+ r'$- \theta_{fr} ^{\circ}C$'
        #if norm == 'frloc':
        var2 = gsw.t_freezing(35,p0,0)
        #for i in np.arange(0,row,1):
        #    #pt_io = np.mean(extract_var(data,'pt_io*_xy',data_type='xy',tval=tval))
        #    var1[int(i),:] = np.subtract(var1[int(i),:],var1[int(i),0])
        varlabel2 = r'$\theta_{fr}(S_\infty) (^{\circ}C)$'
    elif 'ustar' in ops:
        var2,varlable2 = extract_var(data,'us',data_type='ts')
    elif 'geoice' in ops:
        rho = extract_var(data,'rho_ocean',data_type=data_type,
                            tval=tval,xval=xval,yval=yval,zval=zval)
        pt_io = np.mean(extract_var(data,'pt_io*_xy',data_type='xy',tval=tval))
        zmax = np.min(extract_var(data,'z',data_type='pr'))
        pt_far = extract_var(data,'pt',data_type='pr',tval=tval,zval=zmax)
        sa_far = extract_var(data,'sa',data_type='pr',tval=tval,zval=zmax)
        thermal_driving = np.mean(extract_var(data,'thermal_driving',data_type='2d',
                                                tval=tval))
        dT_far = pt_far - pt_io
        drho = beta_S*sa_far - beta_T*dT_far 
        var2 = (g/f) * math.sin(alpha) * (drho/dT_far) * thermal_driving
        varlabel2 = r'$U_{gi}$'
    elif 'geo' in ops:
        # base the geostrophic velocity on the boundary conditions imposed 
        # rather than recalculate
        zmax = np.min(extract_var(data,'z',data_type='pr'))
        u = np.square(extract_var(data,'u',tval=0,zval=zmax))
        v = np.square(extract_var(data,'v',tval=0,zval=zmax))
        var2 = np.sqrt(u,v)
        axis_label = r'$U_g$'
    elif 'far' in ops:
        var2,varlabel2 = extract_var(data,varname,data_type=data_type,
                         tval=tval,zval=[-1000,-1000],ops='tmean',keep='z')
        varlabel2 = varlabel[varlist.index(varname)] + r'_{\infty}'
    if 'norm' in ops:
        var1 = np.divide(var1,var2)
        axis_label = ( varlabel[varlist.index(varname)]+ 
                       '$/$' + varlabel2 )
    elif 'diff' in ops:
        var1 = np.subtract(var1,var2)
        varlabel1 = varlabel1 + r'-'

    axis_label = varlabel1 + varlabel2 + r'\:(' + varunits[varlist.index(varname)] + r')'
    
    return var1,r'$'+axis_label+r'$'

#------------------------------------------------------------------------------
# DERIVED_VAR
#   function to hold derived variable calculations
#   Contains the following varname options
#   'b11'
#   'diss'
#   'dEdt'
#   'eta_star' stability paramter
#   'freezing'
#   'Fbuoy' buoyancy flux
#   'Fshear'
#   'Ftrans'
#   'FbuoyT'
#   'FbuoyS'
#   'haline_driving'
#   'K_surf'
#   'km_eff'
#   'kh_eff'
#   'N2'    buoyancy frequency
#   'prho_gsw'
#   'pt1_t'
#   'rho_gsw'
#   'Rf'    flux richarson number
#   'Ri'    gradient richardson number 
#   'S2'    shear
#   'Sw'
#   'thermal_driving'
#   'u2v2' sum of horizontal variance terms
#   'U'
#   'vel_var_ratio'
#   'wb'
#   'X"' use double quotes at the end to determine deviation using 3d data
#   'z_BL'
#   'z_SL'
#   'zeta'
#   'zmin_w2'
# Inputs:
#   data     loaded NETCDF variable from DATA_1D_PR_NETCDF
#   varname  variable name as defined in this function
#------------------------------------------------------------------------------
def derived_var(data,varname,slice_obj,z_offset = 0.,f = gsw.f(-70.),filedir='',
                data_type='pr',p0=10.,data2=[],data_type2='ts'):
    
    print('derive var: '+varname)
    # data should be profile DATA_1D_PR
    if varname[0] == 'b':
        u2 = data.variables['u*2'][:]
        v2 = data.variables['v*2'][:]
        w2 = data.variables['w*2'][:]
        if varname[1:2] == '11':
            var1 = np.divide(u2,v2) + np.divide(u2,w2) 
        elif varname[1:2] == '22':
            var1 = np.divide(v2,u2) + np.divide(v2,w2) 
        elif varname[1:2] == '33':
            var1 = np.divide(w2,u2) + np.divide(w2,v2) 
        elif varname[1:2] == '12':
            uv = data.variables['u*v*'][:]
            var1 = np.divide(w2,u2) + np.divide(w2,v2) 
    
    elif varname == 'diss':
        var1 = np.add(0,1);
    
    elif varname == 'dEdt':
        var1 = data.variables['e'][:]
        dt = data.variables['time'][1:] - data.variables['time'][0:-1]
        de = data.variables['e'][1:,:] - data.variables['e'][0:-1,:]
        for idx,i in enumerate(dt):
           var1[idx+1,:] = de[idx]/i
    
    elif ( varname[0] == 'd') and (varname[-2] == 'd' and len(varname) > 3):
        num = varname[1:-2]
        den = varname[-1]
        grid = vartype[varlist.index(num)] 
        if (den == 'z') or (den == 'x') or (den == 'y'):
            den,den_type = grid_var(den,data_type=data_type,grid=grid)
        if (den == 't'):
            den == 'time'
            den_type = 't'
        var1,temp = extract_var(data,num,data_type=data_type)
        var2,temp = extract_var(data,den,data_type=den_type)
        if len(var1) == 0:
            return data
        [nt,nden] = np.shape(var1)
        dvar = np.subtract(var1[:,1:],var1[:,:-1])
        dden = np.subtract(var2[1:],var2[:-1])
        var1 = np.multiply(np.ones((nt,nden)),np.nan)
        for i in range(nt):
            var1[i,1:] = np.divide(dvar[i,:],dden)

    elif varname == 'eta_star':
        var1 = ( 1. + np.multiply( xi_N,data.variables['u*'][:]) / 
                      np.multiply( abs(f) * Ri_crit,data.variables['ol'][:]) )**-0.5 
    
    elif varname == 'Fbuoy':
        wsa,temp = extract_var(data,'wsa')
        wpt = data.variables['wpt']
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        var1  = np.multiply(-1*g,
                       ( np.multiply(data.variables['alpha_T'],
                                     wpt,dtype=float)
                        - np.multiply(data.variables['beta_S'], 
                                                  wsa,dtype=float)))
    elif varname == 'FbuoyT':
        wpt = data.variables['wpt'][:]
        var1  = np.multiply(-1*g,
                           ( np.multiply(data.variables['alpha_T'],
                                         wpt,dtype=float)))
    elif varname == 'FbuoyS':
        wsa,temp = extract_var(data,'wsa')
        #wsa = np.divide(data.variables['w*sa*'],data.variables['rho_ocean'])
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        var1  = np.multiply(g,
                            ( np.multiply(data.variables['beta_S'], 
                                          wsa,dtype=float)))
    
    elif varname == 'Ftrans':
        var1 = np.add(data.variables['w*u*u*:dz'],data.variables['w*p*:dz'])
    
    elif varname == 'Fshear':
        dudz,temp = extract_var(data,'dudz')
        dvdz,temp = extract_var(data,'dudz')
        var1 = -1. * ( np.multiply(data.variables['w*u*'],dudz) +
                       np.multiply(data.variables['w*v*'],dvdz) )
    
    #elif varname == 'freezing':
    #    data.variables['freezing'] = np.zeros((np.shape(data.variables['melt*_xy'][:])))
    #    data.variables['freezing'][:] = np.ceil(data.variables['melt*_xy'][:],0.0)

    elif varname == 'haline_driving':
        var1 = np.subtract(data.variables['sa1*_xy'],
                           data.variables['sa_io*_xy'])
    
    elif varname == 'kh_eff':
        var1 = 0

    elif varname == 'km_eff':
        wu = extract_var(data,'wu',data_type=data_type) 
        wv = extract_var(data,'wv',data_type=data_type) 
    
    elif varname == 'K_surf':
        if not ('eta_star' in data.variables.keys()):
            data = derived_var(data,'eta_star')
        var1 = ( 0.4 * 0.0052 * 
               np.square(
                  np.multiply(data.variables['eta_star'][:],
                              data.variables['u*'][:]        )
               )/abs(f) )
    
    elif (varname == 'N2') or (varname == 'N1'):
        rho0,temp = extract_var(data,'rho_ocean',tval = [1,1],zval = [-1000,-1000],
                           data_type=data_type,keep='z') 
        if len(rho0) == 0:
            return data
        drhodz,temp = extract_var(data,'drho_oceandz',data_type=data_type)
        #drhodz[drhodz==0] = np.nan
        #if drhodz.size == 0:
        #    return data
        #shp = np.shape(drhodz) #only applicable for data_type = 'pr'
        #data.variables['N2'] = np.nan * np.ones((nt,nz))
        var1 = ( g / rho0 ) * drhodz 
        #data.variables['N1'] = ( -g / rho0 ) * drhodz 
        #data.variables['N1'] = np.sqrt(( -g / rho0 ) * drhodz)
    
    elif varname == 'prho_gsw':
        pt = data.variables['pt'][:]
        sa = data.variables['sa'][:]
        ct = gsw.CT_from_pt(pt,sa)
        var1 = gsw.sigma0(sa,ct)
        
    elif varname == 'pt1_t':
        data2 = load_data(filedir)
        pt = data2.variables['pt'][:]
        koff = data.variables['k_offset_mcph'][:]
        var1 = np.ones((len(koff),1))
        for i in koff:
            data.variables['pt1_t'][i] = pt[tidx,len(pt) - int(i)]
   
    elif varname == 'Rf':
        Fbuoy,temp = extract_var(data,'Fbuoy',data_type=data_type)
        Fshear,temp = extract_var(data,'Fshear',data_type=data_type)
        var1 = Fbuoy
        var1[:,0] = 0
        var1[:,1:] = np.divide(Fbuoy[:,1:],Fshear[:,1:],dtype=float);
   
    elif varname == 'Ri':
        N2,temp = extract_var(data,'N2',data_type=data_type)
        S2,temp = extract_var(data,'S2',data_type=data_type)
        if (len(N2) == 0) or (len(S2) == 0):
            return data
        S2[S2 == 0] = np.nan;
        #S2 = min(S2,1e-10)
        var1 = np.divide(N2,S2,dtype=float);

    elif varname == 'rho_gsw':
        pt = data.variables['pt'][:]
        sa = data.variables['sa'][:]
        if 'hyp' in data.variables.keys():
            P = extract_var(data,'hyp',data_type='pr')
        else:
            rho0 = 1000 #actual rho0 is needed
            z = extract_var(data,'z',grid=grid_type)
            P = p0 + g*rho0*z # now this is not accurate
        ct = gsw.CT_from_pt(pt,sa)
        if grid_type == '3d' or grid_type == '2d':
            x = extract_var(data,'x',grid=grid_type)
            p,temp = np.meshgrid(P,x)
            var1 = gsw.rho(sa,ct,p)
        else:
            var1 = gsw.rho(sa,ct,P)
    
    elif varname == 'S2':
        dudz,temp = extract_var(data,'dudz',data_type=data_type)
        dvdz,temp = extract_var(data,'dudz',data_type=data_type)
        var1 = ( np.square(dudz,dtype=float) + 
                 np.square(dvdz,dtype=float)   );
    
    elif varname == 'Sw':
        var1 = np.divide(data.variables['w*3'][:],
                         np.sqrt(np.power(data.variables['w*2'],3)))
    
    elif varname == 'thermal_driving':
        if data_type == '2d':
            var1 = np.subtract(data.variables['pt1*_xy'],data.variables['pt_io*_xy'])
    
    elif varname == 'dT':
        if data_type == 'ts':
            data2 = load_data(filedir,data_type='pr')
            pt_io,temp = extract_var(data2,'pt',zval=[1,1],keep='z')
            t,temp = extract_var(data2,'time')
            ts,temp = extract_var(data,'time')
            k,temp = extract_var(data,'k_offset_mcph',teval=[i,i])
            var1 = math.nan*k
            for idx,i in t:
               tval,j = ar.find_nearest(ts,i)
               pt_far,temp = extract_var(data2,'pt',zeval=[z[-1*k],z[-1*k]],keep='z')
               var1[j] = pt_far - pt_io[idx]

    elif varname == 'u2v2':
        var1 = np.add(data.variables['u*2'],
                                        data.variables['v*2'])
                          
    elif varname == 'U':
        u,label = (extract_var(data,'u'))
        v,label = (extract_var(data,'v'))
        var1 = np.sqrt(np.square(u),np.square(v))
    
    elif varname == 'u2v2':
        var1 = np.add(data.variables['u*2'],data.variables['v*2'])
    
    elif varname == 'w"U"0':
        var1 = np.sqrt(np.add(np.square(data.variables['w"u"0']),
                               np.square(data.variables['w"v"0'])))
    elif varname == 'wb':
        wsa,temp = extract_var(data,'wsa')
        wpt = data.variables['wpt']
        #wsa = np.divide(data.variables['wsa'],data.variables['rho_ocean'])
        #wsa = np.divide(wsa,data.variables['rho_ocean'])
        data.variables['wb']  = np.multiply(g,
                                   ( np.multiply(data.variables['alpha_T'],
                                                 wpt,dtype=float)
                                    - np.multiply(data.variables['beta_S'], 
                                                wsa,dtype=float)))
    
    elif varname == 'zeta':
        data.variables['zeta'] = np.ones((len(data.variables['time'][:]),1))
        eta_star = ( 1. + np.multiply( xi_N,data.variables['u*'][:]) / 
                          np.multiply( abs(f) * Ri_crit,data.variables['ol'][:]) )**0.5 
        data.variables['zeta'] = ( abs(f) * z_offset / 
                                   ( np.multiply(eta_star,data.variables['u*'][:]) ) )
    elif varname == 'z_BL_N2max':
        zidx = np.argmax(NN.data, axis=1)
        hb = np.abs(NN.z[zidx].mean())
    
    elif varname == 'zmin_w2':
        zidx = [0,-1]
        for idx,i in enumerate(data.variables['time'][:]):
           w2 = np.min(data.variables['w*2'][idx,zidx[0]:zidx[1]])
    
    elif varname == 'z_SL':
        eta_star = derived_var(data,'eta_star')
        us = extract_var(data,'u*')
        var1 = xi_N * ( np.multiply(eta_star,us)/abs(f))
    
    elif varname == 'z_BL':
        eta_star = derived_var(data,'eta_star')
        us = extract_var(data,'u*')
        var1 = 0.4 * ( np.multiply(eta_star,us))/abs(f)
    elif varname == 'vel_var_ratio':
        u2 = data.variables['u*2'][:]
        v2 = data.variables['v*2'][:]
        w2 = data.variables['w*2'][:]
        var1 = np.divide(2.0*w2,(u2+v2))
    else:
        print(varname,' not available')
    return var1[slice_obj]

def slice_var(data,varname,data_type='pr',tunits='hr',
              tval=[9999,9999],xval=[9999,9999],
              yval=[9999,9999],zval=[9999,9999],grid='sc'):  #,keep=''  
    print('sl: '+varname+', data_type= '+data_type+', grid= '+grid)
    if varname not in data.variables.keys():
        varname = varsvars[varsname.index(varname)][0]
        print('sl_new:'+varname)
    varshape = np.shape(data.variables[varname])
    varaxes = data_dim[datafilelist.index(data_type)].copy()
    
    bool_slice = np.zeros((varshape),dtype=bool)
    
    idx_min = np.zeros(len(varaxes),dtype = int)
    idx_max = np.zeros(len(varaxes),dtype = int)
    for a,axis in enumerate(varaxes):
        if axis == 't':
            if (tval[0] != 9999):
                t,temp = extract_var(data,'time',data_type='ts',ops=tunits);
                [val1,idx_min[a]] = ar.find_nearest(t,tval[0])
                [val1,idx_max[a]] = ar.find_nearest(t,tval[1])
                del t
                if idx_min[a] == idx_max[a]:
                    idx_max[a] = idx_min[a] + 1
            else:
                idx_min[a] = 0
                idx_max[a] = varshape[a]
        elif axis == 'z':
            if (zval[0] != 9999):
                if data_type=='z':
                   z,temp = extract_var(data,varname,data_type=data_type,grid=grid)
                else:
                   z,temp = extract_var(data,'z',data_type=data_type,grid=grid)
                [val1,idx_min[a]] = ar.find_nearest(z,zval[0])
                [val1,idx_max[a]] = ar.find_nearest(z,zval[1])
                if idx_min[a] == idx_max[a]:
                    idx_max[a] = idx_min[a] + 1
            else:
                idx_min[a] = 0
                idx_max[a] = varshape[a]
        elif axis == 'y':
            if (yval[0] != 9999):
                y,temp = extract_var(data,'y',data_type=data_type,grid=grid);
                [val1,idx_min[a]] = ar.find_nearest(y,yval[0])
                [val1,idx_max[a]] = ar.find_nearest(y,yval[1])
                if idx_min[a] == idx_max[a]:
                    idx_max[a] = idx_min[a] + 1
            else:
                idx_min[a] = 0
                idx_max[a] = varshape[a]
        elif axis == 'x':
            if (xval[0] != 9999):
                x,temp = extract_var(data,'x',data_type=data_type,grid=grid) #right grid
                [val1,idx_min[a]] = ar.find_nearest(x,xval[0])
                [val1,idx_max[a]] = ar.find_nearest(x,xval[1])
                if idx_min[a] == idx_max[a]:
                    idx_max[a] = idx_min[a] + 1
            else:
                idx_min[a] = 0
                idx_max[a] = varshape[a]
    #idx_slice = [idx_min,idx_max]
    slice_obj = ()
    for i in range(len(varaxes)):
        slice_obj += (slice(idx_min[i],idx_max[i]),)
    #if (len(varaxes) == 1):
    #    slice_obj = (slice(idx_min[0],idx_max[0]))
    #elif len(varaxes) == 2:
    #    slice_obj = (slice(idx_min[0],idx_max[0]),
    #                 slice(idx_min[1],idx_max[1]))
    #elif len(varaxes) == 3:
    #    slice_obj = (slice(idx_min[0],idx_max[0]),
    #                 slice(idx_min[1],idx_max[1]),
    #                 slice(idx_min[2],idx_max[2]))
    #elif len(varaxes) == 4:
    #    slice_obj[idx_min[0]:idx_max[0],idx_min[1]:idx_max[1],
    #              idx_min[2]:idx_max[2],idx_min[3]:idx_max[3]]=True;
    #del var1
    #shp = np.shape(var);
    #sq = ()
    #for idx,i in enumerate(varaxes):
    #    if (i not in keep) and (shp[idx] == 1):
    #        sq += (int(idx),)
    #if len(sq) > 0:
    #    var = np.squeeze(var,axis=sq);
    #    del sq

    return slice_obj,varaxes

def v_gi(data,f,alpha):
    drho = data.variables['rho_ocean'][:] - data.variables['rho_ocean'][:]
    dT_far = data.variables['pt'][:] - data.variables['pt'][:]
    dT_local = data.variables['pt'][:] - data.variables['pt'][:]
    v_gi = (g/f) * math.sin(alpha) * (drho/dT_far) * dT_local
    return v_gi

def compute_variance(var_in):
    # var_in must have dimensions z,y,x
    shape = np.shape(var_in)
    var_out = var_in.copy()
    for i in range(shape[0]):
        var_mean = np.mean(var[i,:,:])
        var_out[i,:,:] = np.subtact(var_in[i,:,:],var_mean)
    return var_out

