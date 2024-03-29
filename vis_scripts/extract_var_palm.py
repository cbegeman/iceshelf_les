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
from var_param_palm import *
from run_table_mod import value_from_namelist

def load_data(filedir,data_type='pr',av=False,coupled=False):
    if data_type == 'parameter':
        data_type = 'ts'
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
    #print(f'grid variable is {var} for variable type {data_type}')
    return var,data_type

def extract_var(data,var_name,data_type='pr',ops=[],
                grid='',keep='t',slice_obj=[],
                tval=[-9999.,-9999.],xval=[9999.,9999.],
                yval=[9999.,9999.],zval=[9999.,9999.],
                tunits = 'hr', tav = 0., 
                varaxes = ['t'], p0 = 10., data_dir=''):
    print('extract_var',var_name,data_type,slice_obj)
    if data_type == 'pr':
        tav_native = 1.
    else:
        tav_native = 0.
    if tav < tav_native:
       tav = 0.
    if data_type == 'parameter':
        #slice_obj_der,_ = slice_var(data,source_var[varlist.index(var_name)],data_type='ts')
        var1 = derived_var(data,var_name,[],filedir=data_dir,p0=p0,
                           data_type=data_type)
        axis_label = varlabel[varlist.index(var_name)]
        if (varunits[varlist.index(var_name)] != r''):
           axis_label = axis_label + r'\:(' + varunits[varlist.index(var_name)] + r')'
        return var1,r'$'+axis_label+r'$'

    if (var_name == 'z') or (var_name == 'x') or (var_name == 'y'):
        grid = vartype[varlist.index(var_name)] 
        var_name,data_type = grid_var(var_name,data_type=data_type,grid=grid)
    if (var_name == 'zE'):
        grid = 'z'#vartype[varlist.index(var_name)] 
    if slice_obj == [] and data_type != 'parameter':
        print('____slice_var called in extract_var for ',var_name)
        slice_obj,varaxes = slice_var(data,source_var[varlist.index(var_name)],data_type=data_type,
                                      tval = tval,
                                      xval=xval, yval=yval, zval=zval, grid=grid)
    if var_name in dervar:
        print('____derived_var called in extract_var for ',var_name)
        slice_obj_der,_ = slice_var(data,source_var[varlist.index(var_name)],data_type=data_type)
        var1 = derived_var(data,var_name,
                           slice_obj_der,
                           filedir=data_dir,p0=p0,
                           data_type=data_type,tlim=[tval[0]-tav/2.,tval[1]+tav/2.])
        #if np.shape(data_type)==np.shape(slice_obj):
    else:
        if not(var_name in data.variables.keys()):
            print('variable ' + var_name + ' not in file')
            return [],''
        var1 = data.variables[var_name][:]
    
    if tav > 0. and 't' in varaxes:
        print('____time average in extract_var')
        #slice_obj_tav,_ = slice_var(data,source_var[varlist.index(var_name)],data_type=data_type,
        #                              tval = [tval[0] - tav/2., tval[1] + tav/2.],
        #                              xval=xval, yval=yval, zval=zval, grid=grid)
        t = np.divide(data['time'][:],3600.)
        #t,_ = extract_var(data,'time')#,slice_obj=slice_obj[varaxes.index('t')])
        ti = tval[0]
        tmask = np.multiply(t>(ti - tav/2.),t<=(ti + tav/2.))
        #print('____extract_var time mean: '+var_name+' shp = ',np.shape(var1),'tmask',t[tmask])
        var1_shape = np.shape(var1)
        var1 = var1.data
        if len(var1_shape) == 2:#data_type == 'pr':
            #print('____time average in extract_var for pr')
            #print('shape tmask',np.shape(tmask))
            #print('shape var1 tav',np.shape(var1))
            #print('tmask',t[tmask])
            #print('type var1',type(var1))
            #print(np.shape(var1[tmask,:]),varaxes.index('t'))
            _,nz = np.shape(var1)
            for i in range(nz):
                #var1[np.argmin(np.abs(t-ti)),i] = np.mean(np.mean(var1[tmask,i]))
                var1[:,i] = np.nanmean(var1[tmask,i])
                #var1[np.argmin(np.abs(t-ti)),i] = np.mean(var1[tmask,i])
        elif len(var1_shape) == 1:
            if data_type == 'ts':
                var1[:] = np.divide(np.sum(np.multiply(var1[tmask],data['dt'][tmask])),np.sum(data['dt'][tmask]))
            else:
                #var1[np.argmin(np.abs(t-ti))] = np.mean(var1[tmask])
                var1[:] = np.nanmean(var1)
        else: 
            print('Unsupported data type in extract_var with time averaging')
            return
        #temp = var1
        #tmask = np.multiply(t<(tval[1] - tav/2.), t>=(tval[0] + tav/2.))
        #tmask2 = np.zeros(np.shape(var1),dtype=bool)
        #if ~np.any(tmask):
        #    tmask[np.argmin(np.abs(t-(tval[0]+tav/2.)))] = True
        #for i,ti in enumerate(t):
        #    if tmask[i]:
        #        tmask2 = np.multiply(t>(ti - tav/2.),t<=(ti + tav/2.))
        #        var1[i,:] = np.mean(temp[tmask2,:],axis=varaxes.index('t'))
        #slice_obj,varaxes = slice_var(data,source_var[varlist.index(var_name)],data_type=data_type,
        #                              tval=tval,xval=xval,yval=yval,zval=zval,grid=grid)
    var1 = var1[slice_obj]
    #print('extract_var: '+var_name+', data_type='+data_type+', grid=',grid,'slice_obj=',slice_obj,'shape=',np.shape(var1))
    
    shp = np.shape(var1);
    sq = ()
    for idx,i in enumerate(varaxes):
        if (i not in keep) and (shp[idx] == 1):
            sq += (int(idx),)
    if len(sq) > 0:
        var1 = np.squeeze(var1,axis=sq);
        del sq
    
    # perform operations on var 
    varlabel1 = varlabel[varlist.index(var_name)]
    varlabel2 = r''
    if var_name[0:2] == 'pt':
        var1 = np.subtract(var1,K0)
        varunits[varlist.index(var_name)] = r'^{\circ}C'
    elif 'melt' in var_name and tunits=='hr':
        var1 = np.multiply(var1,s_yr)
        varunits[varlist.index(var_name)] = r'm\:yr^{-1}'
    if 'shf' in ops:
        var2,varlabel2 = extract_var(data,'shf',data_type='ts')
    elif 'rho0' in ops:
        var2,varlabel2 = extract_var(data,'rho_ocean',
                         tval=[0,0],zval=[0,0],data_type='pr')
    elif ('mean' in ops) and len(var1) > 1:
        # issue with using axis formulation for dimension 1 arrays
        var2 = np.mean(var1);#,axis = var_axes.index(ops[0]))
        varlabel2 = r'\overline{'+varlabel[varlist.index(var_name)]+r'}'
        if ('diff' not in ops) and ('norm' not in ops):
           var1 = var2.copy() # if no further operation, output mean
    elif 'std' in ops:
        var2 = np.std(var1,axis=varaxes.index(ops[0]))
        varlabel2 = r'\sigma_{'+varlabel[varlist.index(var_name)]+'}'
        if ('diff' not in ops) and ('norm' not in ops):
           var1 = var2.copy()
    elif 'sasws' in ops:
        var2,varlabel2 = extract_var(data,'sasws',data_type='ts')
    elif 'ws3' in ops:
        var2,varlabel2 = extract_var(data,'NORM_ws3',tval=tval,data_type='pr')
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
        thermal_driving = np.mean(extract_var(data,'thermal_driving_i',data_type='2d',
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
        var2,varlabel2 = extract_var(data,var_name,data_type=data_type,
                         tval=tval,zval=[-1000,-1000],ops='tmean',keep='z')
        varlabel2 = varlabel[varlist.index(var_name)] + r'_{\infty}'
    if 'norm' in ops:
        var1 = np.divide(var1,var2)
        axis_label = ( varlabel[varlist.index(var_name)]+ 
                       '$/$' + varlabel2 )
    elif 'diff' in ops:
        var1 = np.subtract(var1,var2)
        varlabel1 = varlabel1 + r'-'

    if tunits=='hr' and var_name == 'time':
        var1 = np.divide(var1,3600.)
        varunits[varlist.index(var_name)] = 'h'
    
    axis_label = varlabel1 + varlabel2
    if (varunits[varlist.index(var_name)] != r''):
       axis_label = axis_label + r'\:(' + varunits[varlist.index(var_name)] + r')'
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
def derived_var(data,varname,slice_obj_input,z_offset = 0.,f = gsw.f(-70.),filedir='',
                tlim=[-9999.,-9999.], data_type='pr',p0=10.,data2=[],data_type2='ts'):
    #print('derive_var: '+varname+', data_type= '+data_type,'slice_obj=',slice_obj_input)

    # data should be profile DATA_1D_PR
    if varname == 'sin_alpha':
        base_dir,_ = filedir.split('RUN')
        alpha = value_from_namelist(base_dir,'alpha_surface')[0]
        var1 = [math.sin(math.pi*alpha/180)]
    elif varname[0] == 'b':
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
        dedt,_ = extract_var(data,'dEdt',data_type=data_type,data_dir=filedir)
        Fshear,_ = extract_var(data,'Fshear',data_type=data_type,data_dir=filedir)
        Ftrans,_ = extract_var(data,'Ftrans',data_type=data_type,data_dir=filedir)
        Fbuoy,_ = extract_var(data,'Fbuoy',data_type=data_type,data_dir=filedir)
        var1 = np.subtract(np.add(Fshear,Ftrans,Fbuoy),dedt)
    
    elif varname == 'zE':
        t = data.variables['time'][:]/3600.
        ti = np.argmin(np.abs(t-50.))
        #var1,_ = extract_var(data,'z',slice_obj=slice_obj)#,grid='zu'
        K_ref,_ = extract_var(data,'km_eff',data_type=data_type)#,slice_obj=slice_obj,tav=tav)
        K_ref = K_ref.data
        print('shp(K_ref)',np.shape(K_ref))
        var1 = np.sqrt(2*np.nanmean(K_ref[ti-12:ti+1,-80:]))/np.sqrt(abs(f))
        #var1 = np.divide(var1,dE)

    elif varname == 'dEdt':
        var1 = data.variables['e*'][:] 
        dt = data.variables['time'][1:] - data.variables['time'][:-1]
        de = data.variables['e*'][1:,:] - data.variables['e*'][:-1,:]
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
        dden = np.maximum(1e-10,dden)
        var1 = np.multiply(np.ones((nt,nden)),np.nan)
        for i in range(nt):
            var1[i,1:] = np.divide(dvar[i,:],dden)

    elif varname == 'eta_star':
        var1 = ( 1. + np.multiply( xi_N,data.variables['u*'][:]) / 
                      np.multiply( abs(f) * Ri_crit,data.variables['ol'][:]) )**-0.5 
    
    elif varname == 'gamma_T_2m':
        #print('____gamma_T_2m slice_obj in derived_var',slice_obj_input)
        #gammaT = derived_var(data,'gammaT',slice_obj,data_type=data_type)
        data2 = load_data(filedir,data_type='ts')
        time,_ = extract_var(data2,'time',data_type='ts')#,slice_obj=slice_obj_input)
        dt,_ = extract_var(data2,'dt',data_type='ts')#,slice_obj=slice_obj_input)
        melt_t,_ = extract_var(data2,'melt',data_type='ts',tunits='s')#,slice_obj=slice_obj_input)
        pt0_t,_ = extract_var(data2,'pt(0)',data_type='ts')#,slice_obj=slice_obj_input)

        t0,_ = extract_var(data,'time',data_type='pr',slice_obj=slice_obj_input[0])
        t_pr,_ = extract_var(data,'time',data_type='pr',tval=tlim)
        pt,_ = extract_var(data,'pt',data_type='pr',zval=[-2.,-2.],tval=tlim)
        u,_ = extract_var(data,'U',data_type='pr',zval=[-2.,-2.],tval=tlim)
        us = (0.003)**(0.5)*u
        
        var1 = np.zeros(np.shape(t0))
        melt = np.zeros(np.shape(t_pr))
        dpt = np.zeros(np.shape(t_pr))
        for i,t in enumerate(t_pr[:-1]):
            tmask = np.multiply(time<t_pr[i+1],time >= t_pr[i])
            melt[i] = np.divide(np.sum(np.multiply(melt_t[tmask],dt[tmask])),np.sum(dt[tmask]))
            dpt[i] = pt[i] - np.divide(np.sum(np.multiply(pt0_t[tmask],dt[tmask])),np.sum(dt[tmask]))
        tmask = np.multiply(time<t_pr[-1]+1,time >= t_pr[-1])
        melt[-1] = np.mean(melt_t[tmask])
        dpt[-1] = pt[-1] - np.mean(pt0_t[tmask])
        var1[:] = (L_i/cpw) * np.divide(np.mean(melt),np.multiply(np.mean(us),np.mean(dpt)))
    
    elif varname == 'gamma_T':
        #gammaT = derived_var(data,'gammaT',slice_obj,data_type=data_type)
        melt,_ = extract_var(data,'melt',data_type=data_type,tunits='s')
        pt0,_ = extract_var(data,'pt(0)',data_type=data_type)
        pt_zmo,_ = extract_var(data,'pt(z_mo)',data_type=data_type)
        us,_ = extract_var(data,'u*',data_type=data_type)
        var1 = (L_i/cpw) * np.divide(melt, 
                                     ( np.multiply(us,
                                                   np.subtract(pt_zmo,pt0) ) ) )
    elif varname == 'gammaT':
        melt,_ = extract_var(data,'melt',data_type=data_type,tunits = 's')
        pt0,_ = extract_var(data,'pt(0)',data_type=data_type)
        pt_zmo,_ = extract_var(data,'pt(z_mo)',data_type=data_type)
        var1 = (L_i/cpw) * np.divide(melt, 
                          ( np.subtract(pt_zmo,pt0) ) )
    elif varname == 'Fbuoy_u':
        usa = np.add(data.variables['u*sa*'],data.variables['u"sa"'])
        upt = np.add(data.variables['u*pt*'],data.variables['u"pt"'])
        base_dir,_ = filedir.split('RUN')
        alpha = value_from_namelist(base_dir,'alpha_surface')[0]
        var1  = np.multiply(g*math.sin(math.pi*alpha/180),
                       ( np.multiply(data.variables['alpha_T'],
                                     upt,dtype=float)
                        - np.multiply(data.variables['beta_S'], 
                                                  usa,dtype=float)))
    elif varname == 'Fbuoy_w':
        wsa = data.variables['wsa']
        wpt = data.variables['wpt']
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        var1  = np.multiply(g,
                       ( np.multiply(data.variables['alpha_T'],
                                     wpt,dtype=float)
                        - np.multiply(data.variables['beta_S'], 
                                                  wsa,dtype=float)))
    elif varname == 'Fbuoy':
        usa = np.add(data.variables['u*sa*'],data.variables['u"sa"'])
        upt = np.add(data.variables['u*pt*'],data.variables['u"pt"'])
        base_dir,_ = filedir.split('RUN')
        alpha = value_from_namelist(base_dir,'alpha_surface')[0]
        fu = np.multiply(g*math.sin(math.pi*alpha/180),
                       ( np.multiply(data.variables['alpha_T'],
                                     upt,dtype=float)
                        - np.multiply(data.variables['beta_S'], 
                                                  usa,dtype=float)))
        wsa = data.variables['wsa']
        wpt = data.variables['wpt']
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        wsa = np.divide(wsa,data.variables['rho_ocean'])
        fw  = np.multiply(g,
                       ( np.multiply(data.variables['alpha_T'],
                                     wpt,dtype=float)
                        - np.multiply(data.variables['beta_S'], 
                                                  wsa,dtype=float)))
        var1 = np.add(fw,fu)

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
        wu,_ = extract_var(data,'wu',data_type=data_type) 
        wv,_ = extract_var(data,'wv',data_type=data_type) 
        dudz,_ = extract_var(data,'dudz',data_type=data_type)
        dvdz,_ = extract_var(data,'dvdz',data_type=data_type)
        wU = np.sqrt(np.add(np.square(wu),np.square(wv)))
        dUdz = np.sqrt(np.add(np.square(dudz),np.square(dvdz)))
        var1 = np.zeros(np.shape(wu))
        var1 = np.divide(wU,dUdz)
    
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
        Fbuoy,temp = extract_var(data,'Fbuoy',data_type=data_type)#, tav=tav)
        Fshear,temp = extract_var(data,'Fshear',data_type=data_type)#, tav=tav)
        var1 = Fbuoy
        var1[:,0] = 0
        var1[:,1:] = np.divide(Fbuoy[:,1:],Fshear[:,1:],dtype=float);
   
    elif varname == 'Ri':
        N2,temp = extract_var(data,'N2',data_type=data_type)#, tav=tav)
        S2,temp = extract_var(data,'S2',data_type=data_type)#, tav=tav)
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
        dudz,temp = extract_var(data,'dudz',data_type=data_type)#, tav=tav)
        dvdz,temp = extract_var(data,'dudz',data_type=data_type)#, tav=tav)
        var1 = ( np.square(dudz,dtype=float) + 
                 np.square(dvdz,dtype=float)   )
    
    elif varname == 'Sw':
        var1 = np.divide(data.variables['w*3'][:],
                         np.sqrt(np.power(data.variables['w*2'],3)))
    
    elif varname == 'thermal_driving_i':
        if data_type == '2d':
            var1 = np.subtract(data.variables['pt1*_xy'],data.variables['pt_io*_xy'])
        if data_type == 'ts':
            pt0,_ = extract_var(data,'pt(0)',data_type=data_type)
            pt_zmo,_ = extract_var(data,'pt(z_mo)',data_type=data_type)
            var1 = np.subtract(pt_zmo,pt0)
    
    elif varname == 'thermal_driving_infty':
        base_dir,_ = filedir.split('RUN')
        tf = gsw.t_freezing(35,800,0)
        pt_freezing = gsw.pt0_from_t(35,tf,800)
        pt_surf = value_from_namelist(
                       base_dir,
                       'pt_surface')[0] - K0 - pt_offset
        var1 = [pt_surf - pt_freezing]
        #pt_far = extract_var(data,'pt',data_type='pr',tval=tval,zval=zmax)
        #pt_far = value_from_namelist(
        #           '/lustre/scratch3/turquoise/cbegeman/palm/jobs/'+filedir,
        #           'pt_surface')[0] - K0
        #if data_type == '2d':
        #    var1 = np.subtract(data.variables['pt_io*_xy'],pt_far)
        #if data_type == 'ts':
        #    pt0,_ = extract_var(data,'pt(0)',data_type=data_type)
        #    var1 = np.subtract(pt0,pt_far)
    
    elif varname == 'Umax':
        U,_ = extract_var(data,'U',zeval=[0,9999])#,tval=[t,t])
        var1 = np.max(U)

    elif varname == 'Umax_i':
        for t in range(tval[0],tval[0]+tav):
            U,_ = extract_var(data,'U',zval=[0,9999])#,tval=[t,t])
            var1 = max(var1,U)
    elif varname == 'u*':
        u,_ = extract_var(data,'u',data_type=data_type)#,zval=[0,9999],tval=[t,t])
        var1 = np.subtract(u,np.mean(u))
    elif varname == 'v*':
        u,_ = extract_var(data,'v',data_type=data_type)#,zval=[0,9999])#,tval=[t,t])
        var1 = np.subtract(u,np.mean(u))
    elif varname == 'U*':
        u,_ = extract_var(data,'u',data_type=data_type)#,zval=[0,9999])#,tval=[t,t])
        u_prime = np.subtract(u,np.mean(u))
        v,_ = extract_var(data,'v',data_type=data_type)#,zval=[0,9999])#,tval=[t,t])
        v_prime = np.subtract(v,np.mean(v))
        var1 = np.sqrt(np.add(np.sq(u_prime),np.sq(v_prime)))

    elif varname == 'dT':
        if data_type == 'ts':
            data2 = load_data(filedir,data_type='pr')
            pt_io,temp = extract_var(data2,'pt',zval=[1,1],keep='z')
            t,temp = extract_var(data2,'time')
            ts,temp = extract_var(data,'time')
            k,temp = extract_var(data,'k_offset_mcph',tval=[i,i])
            var1 = math.nan*k
            for idx,i in t:
               tval,j = ar.find_nearest(ts,i)
               pt_far,temp = extract_var(data2,'pt',zeval=[z[-1*k],z[-1*k]],keep='z')
               var1[j] = pt_far - pt_io[idx]

    elif varname == 'u2v2':
        var1 = np.add(data.variables['u*2'],
                                        data.variables['v*2'])
                          
    elif varname == 'U':
        u,label = extract_var(data,'u')#,slice_obj=slice_obj_input)
        v,label = extract_var(data,'v')#,slice_obj=slice_obj_input)
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
    elif varname == 'w*2_max':
        w2 = data.variables['w*2'][:]
        var1 = np.amax(w2,axis=1)
    
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
        var1 = np.divide(2.0*w2,np.maximum(w2/5.,(u2+v2)))
    else:
        print(varname,' not available')
    return var1

def slice_var(data,varname,data_type='pr',tunits='hr',
              tval=[-9999,-9999],xval=[9999,9999],
              yval=[9999,9999],zval=[9999,9999],grid='sc'): 
    #print(f'slice {varname}')
    if varname in dervar:
        varname = source_var[varlist.index(varname)]
        #print('slice_var source_var: ',varname)
    elif varname not in data.variables.keys():
        varname = varsvars[varsname.index(varname)][0]
    varshape = np.shape(data.variables[varname])
    varaxes = data_dim[datafilelist.index(data_type)].copy()
    if varname == 'time':
        varaxes = ['t']
    #print(f'slice_var: {varname},{varshape},{varaxes},{data_type}')
    
    bool_slice = np.zeros((varshape),dtype=bool)
    idx_min = np.zeros(len(varaxes),dtype = int)
    idx_max = np.zeros(len(varaxes),dtype = int)
    for a,axis in enumerate(varaxes):
        if axis == 't':
            if (tval[0] != -9999):
                t,temp = extract_var(data,'time',data_type=data_type,tunits=tunits);
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
    #print(f'{varname} sliced with {varaxes}')
    #slice_obj = [slice(idx_min[i],idx_max[i]) for i in range(len(varaxes))]
    for i in range(len(varaxes)):
        slice_obj += (slice(idx_min[i],idx_max[i]),)
        #print(f'Indexing {varaxes[i]} from {idx_min[i]} to {idx_max[i]}')
    #if (len(varaxes) == 1):
    #    slice_obj = slice(idx_min[0],idx_max[0])
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

