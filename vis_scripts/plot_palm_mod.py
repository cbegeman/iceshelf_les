#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 07:52:34 2019

@author: cbegeman
"""
import math,os
#from scipy.optimize import curve_fit
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib as pltlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as tick
import my_array_mod as ar
import my_string_mod as strop
import convert_NCAR_mod as ncar
from extract_var_palm import *
from matplotlib import rc
from run_table_mod import value_from_namelist
import cmocean
import var_param_palm as pv
from plot_param_palm import *

pltlib.rc_file('rcparams.txt', use_default_template=True)
#rc('text', usetex=True)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('xtick',labelsize=14)
#plt.rc('ytick',labelsize=14)
#plt.rc('axes.formatter',useoffset=False)
        #ax.xaxis.get_major_formatter().set_powerlimits((pw_min,pw_max))
        #ax.yaxis.get_major_formatter().set_powerlimits((pw_min,pw_max))

#------------------------------------------------------------------------------
# END_TIME
# Inputs:
#   runfile path to RUN_CONTROL
# Outputs: 
#   end time of run in seconds of model time
#------------------------------------------------------------------------------
def end_time(runfile): 
   runcontrol = open(runfile,'r')
   linelist = runcontrol.readlines()
   hr = int(linelist[len(linelist)-1][11:13]) 
   minute = int(linelist[len(linelist)-1][14:16]) 
   sec = int(linelist[len(linelist)-1][17:19])
   run_end = hr * 3600. + minute * 60. + sec
   return run_end

#------------------------------------------------------------------------------
# FETCH_RUN_NAMES
# Inputs:
#   base_run path to run directory containing sample namelist
#   sens_var variable that is allowed to vary between runs
# Outputs: 
#   cell array of run directories that are similar to base_run except they have
#   a different value for sens_var
#------------------------------------------------------------------------------
def fetch_run_names(base_run,sens_var,namelist = 'test_oceanml_p3d'):
   idx = [pos for pos, char in enumerate(base_run) if char == '/']
   base_dir = base_run[:idx[-3]+1]
   run_dir = base_run[idx[-3]+1:idx[-2]+1]
   f = open(base_dir + run_dir + namelist,'r')
   dirs = os.listdir(base_dir)
   for diri in dirs:
      if ( (diri.startswith(base_dir)) and 
           (diri != run_dir) and
           os.path.isfile(base_dir + diri + '/' + namelist) ) :
         f2 = open(base_dir + diri + '/' + namelist,'r')
         fdiff = difflib.unified_diff(f,f2)
         #print((list(fdiff)))#if (len(fdiff) > 0) and (len(fdiff) <= 2):
         #   print(diri)
         for line in fdiff:
               print(line)
         f2.close()
   f.close()
   run_names = []
   return run_names,sens_val

#------------------------------------------------------------------------------
# MELT_STATS
#   returns the cumulative melt over the simulation period (default) or a user-
#   defined time period
# Inputs:
#   filedir  path to run directory
#   tlim     vector length 2, start time and end time in hours
# Outputs: 
#   cumulative melt in m
#------------------------------------------------------------------------------
def var_stats(filedir, varname, tlim = [9999.,9999.]):
    data1 = load_data(filedir,data_type='ts')
    t,t_label = extract_var(data1,'time',tunits='hr',data_type='ts',tval=tlim)
    if len(t) < 1:
        return -9999.,0.
    dt = t
    dt[1:] = t[:-1] - t[1:]
    dt[0] = dt[1]
    var1,temp = extract_var(data1,varname,tunits='hr',data_type='ts',tval=tlim)
    mean = np.sum(np.multiply(var1,dt))/(np.sum(dt))
    sd = np.var(var1)
    return mean,sd

#------------------------------------------------------------------------------
# PLOT_2D_XY
#   creates a plot of PALM variable field from XY plane output data
# Inputs:
#   filedir  path to run directory
#   runname  shorthand label for run
#   plotvar  cell array of variable names to be plotted 
#   teval    vector of times to plot field in hours (optional). Default is the 
#            first and last model output times. 
#   zeval    model depth at which to output (optional)
#            TODO check whether this is stored in DATA_2D_XY_NETCDF
#   clim     axis limits for variable in plotvar (optional). 
#            Currently the same limits are applied to every variable in plotvar
#   clim_alltime logical (optional). if true, defines axis limits for variable  
#            in plotvar at the maximum and minimum value of that variable for 
#            all model output times
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#------------------------------------------------------------------------------
def plot_2d_xy(filedir, runname, plotvar, teval = [9999.], zeval = [0],
            tunits = 'hr', clim = [-9999.,-9999.], clim_alltime = False, 
            outputdir = [], printformat = 'png', av=False,prc=5):

    for j in plotvar:

        data1 = load_data(filedir,data_type='2d',av=av)
        x,x_axis_label = extract_var(data1,'x',data_type = '2d',
                                        grid=pv.vartype[pv.varlist.index(j)])
        y,y_axis_label = extract_var(data1,'y',data_type = '2d',
                                        grid=pv.vartype[pv.varlist.index(j)])
        x_mesh, y_mesh = np.meshgrid(x,y)
        
        z,z_axis_label = extract_var(data1,'z',data_type = '2d',keep='z',
                                        grid=pv.vartype[pv.varlist.index(j)])
        [zval,zidx] = ar.find_nearest(z,zeval)
                    
        if clim_alltime:
            var_allt,temp = extract_var(data1,j,data_type = '2d')
            cmin = np.min(var_allt)
            cmax = np.max(var_allt)
            del var_allt

        for i in teval:
            t,t_label = extract_var(data1,'time',ops=tunits,data_type='ts',
                                       keep='t',tval=[i,i])
            print(t) 
            var1,cbar_label = extract_var(data1,j,data_type = '2d',tval=[i,i])
            fig = plt.figure()
            plt.set_cmap(pv.varcmap[pv.varlist.index(j)])
            plt.pcolor(x_mesh,y_mesh,var1)
            plt.xlabel(x_axis_label,fontsize = fs)
            plt.ylabel(y_axis_label,fontsize = fs)
                    
            cbar = plt.colorbar()
            cbar.set_label(cbar_label)
            climits = [0.,0.]
            if clim[0] == -9999.:
                if clim_alltime:
                    climits[0] = np.min(var_allt)
                    climits[1] = np.max(var_allt)
                elif prc>0:
                    climits[0] = np.percentile(var1,prc) 
                    climits[1] = np.percentile(var1,100-prc) 
                else:
                    climits[0] = np.min(var1)
                    climits[1] = np.max(var1)
                if j == 'u' or j == 'v' or j == 'w' or j == 'melt*_xy':
                    climits[0] = -1.*max(abs(climits[0]),abs(climits[1]))
                    climits[1] =     max(abs(climits[0]),abs(climits[1]))
            else:
                climits = clim
            plt.clim(climits)
            
            cbar.set_label(cbar_label,fontsize = fs)
            plt.title(runname + ' ' + t_label + ' = ' + str(int(t)),fontsize = fs)
                    
            if len(outputdir) == 0:
                outputdir = filedir
            name = pv.varname[pv.varlist.index(j)] + '_'
            filenamesave = ( name + 'xy_z' + str(int(zval)) + 
                             '_t' + str(int(t)) + '.' + printformat )
            plt.savefig(outputdir+filenamesave)
            plt.close()
            print(outputdir + filenamesave)

            fig = plt.figure()
            counts,bin_edges = np.histogram(var1,bins='auto')
            cdf = np.cumsum(counts)
            plt.plot(bin_edges[1:],cdf)
            name = pv.varname[pv.varlist.index(j)] + '_'
            filenamesave = ( name + 'xy_z' + str(int(zval)) + 
                             '_t' + str(int(t)) + '_hist.' + printformat )
            plt.savefig(outputdir+filenamesave)
            plt.close()
             
        data1.close()

#------------------------------------------------------------------------------
# PLOT_XY
#   creates a plot of PALM variable field from XY from 3D output data at user-
#   specified depth
# Inputs:
#   filedir  path to run directory
#   runname  shorthand label for run
#   plotvar  cell array of variable names to be plotted 
#   teval    vector of times to plot field in hours (optional). Default is the 
#            first and last model output times. 
#   zeval    model depth at which to output (optional)
#   clim     axis limits for variable in plotvar (optional). 
#            Currently the same limits are applied to every variable in plotvar
#   clim_alltime logical (optional). if true, defines axis limits for variable  
#            in plotvar at the maximum and minimum value of that variable for 
#            all model output times
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#------------------------------------------------------------------------------
def plot_3d_slice(filedir, runname, plotvar, runlabel = '', 
                  teval = [9999.], 
                  ops='', tunits = 'hr', slice_dim = '',
                  xval = [9999.,9999.], yval = [9999.,9999.], zval = [9999.,9999.],
                  clim = [-9999.,-9999.], clim_alltime = False, 
                  outputdir = [], printformat = 'png', av=False,prc = 1, figsize=(6.4,4.8),
                  filecmp = '', runnamecmp = '',overwrite=True):
    
    if len(plotvar)<1:
        print('No variables listed to be plotted')
        return
    data1 = load_data(filedir,data_type='3d',av=av)
    if runlabel == '':
        runlabel = runname
    if filecmp != '':
        runcmp = True
        datacmp = load_data(filecmp,data_type='3d',av=av)
    else:
        runcmp = False

            
    if zval[0] != 9999.:
       slice_dim = 'z'
    if yval[0] != 9999.:
       slice_dim = 'y'
    if xval[0] != 9999.:
       slice_dim = 'x'

    x, x_axis_label = extract_var(data1,'x',data_type = '3d',xval=xval,
                                  grid=pv.vartype[pv.varlist.index(plotvar[0])],keep='x')
    y, y_axis_label = extract_var(data1,'y',data_type = '3d',yval=yval,keep='y',
                                  grid=pv.vartype[pv.varlist.index(plotvar[0])])
    z, z_axis_label = extract_var(data1,'z',data_type = '3d',zval=zval,keep='z',
                                  grid=pv.vartype[pv.varlist.index(plotvar[0])])
    
    if slice_dim == 'z':
       x_mesh, y_mesh = np.meshgrid(x,y)
       slval = z#[zidx]
    if slice_dim == 'y':
       x_mesh, y_mesh = np.meshgrid(x,z)
       y_axis_label = z_axis_label
       slval = x
    if slice_dim == 'x':
       x_mesh, y_mesh = np.meshgrid(y,z)
       x_axis_label = y_axis_label
       y_axis_label = z_axis_label
       slval = x
    
    for idx,j in enumerate(plotvar):
        
        for i in teval:
            
            tval,t_label = extract_var(data1,'time',data_type ='ts',
                                       keep='t',tunits=tunits,tval=[i,i])
            
            if len(outputdir) == 0:
                outputdir = filedir
            filenamesave = strop.rm_palm_char(j)
            if runcmp:
                filenamesave = filenamesave + '_' + runnamecmp
            if slice_dim == 'z':
                filenamesave = filenamesave + 'xy_z' + str(int(abs(zval[0]))) 
            elif slice_dim == 'y':
                filenamesave = filenamesave + 'xz_y' + str(int(yval[0]))
            if slice_dim == 'x':
                filenamesave = filenamesave + 'yz_x' + str(int(xval[0]))
            if zval[0] != 9999.:
                filenamesave = filenamesave + '_zmax' + str(int(abs(zval[0]))) 
            filenamesave = filenamesave + '_t' + str(int(tval)) + '.' + printformat 
            if os.path.exists(outputdir + filenamesave) and not overwrite:
                print('skipping existing file: '+ outputdir + filenamesave)
                continue 
    
            if runcmp:
                tcmp,no_label = extract_var(datacmp,'time',tunits=tunits,tval=[i,i])
            
            var1, cbar_label = extract_var(data1,j,data_type = '3d',tval=[i,i],
                                           tunits=tunits,xval=xval,yval=yval,zval=zval)
            #print(xval,yval,zval)
            #print(np.shape(x_mesh))
            #print(np.shape(y_mesh))
            #print(np.shape(var1))
            if(len(np.shape(var1))==3):
                var1 = var1[0,:,:]
            #if j == 'u' or j == 'v':
            #    var1 = np.subtract(var1,np.mean(var1))
            #    cbar_label = r'$v^{\prime} \: (m\:s^{-1})$'
            if runcmp:
                [tval2,tidx2] = ar.find_nearest(tcmp,i)
                varcmp, no_label = extract_var(datacmp,j,data_type = '3d',tval=[i,i],
                                               grid=pv.vartype[pv.varlist.index(j)]);
                                                 #tidx=tidx2,xidx=xidx,yidx=yidx,zidx=zidx)
                # check that grids have the same resolution
                varshp = np.shape(var1)
                varshp2 = np.shape(varcmp)
                if (varshp != varshp2):
                    print('Datasets have different sizes.')
                    return
                    #TODO add interpolation
            
                var1 = np.subtract(var1,varcmp)

            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
            if pv.varscale[pv.varlist.index(j)]=='log':
                if np.size(var1[var1>0])>0:
                    logmin = np.min(var1[var1>0])
                    var1 = np.maximum(var1,logmin)
                plt.pcolor(x_mesh,y_mesh,var1,
                           norm=colors.LogNorm(),
                           cmap=pv.varcmap[pv.varlist.index(j)])
            else:
                plt.pcolor(x_mesh,y_mesh,var1,
                           cmap=pv.varcmap[pv.varlist.index(j)])
            plt.xlabel(x_axis_label,fontsize = fs)
            plt.ylabel(y_axis_label,fontsize = fs)
            ax.set_aspect('equal', adjustable='box')
            if slice_dim != 'z':
               ax.set_yticks([0,-20,-40])

            cbar = plt.colorbar()
            cbar.set_label(cbar_label,fontsize = fs)
            climits = [0.,0.]
            if clim[0] != -9999.:
                climits = clim.copy()
            else:
                if clim_alltime:
                    climits[0] = np.min(var_allt)
                    climits[1] = np.max(var_allt)
                else:
                    #climits[0] = np.percentile(var1,prc) 
                    #climits[1] = np.percentile(var1,100-prc) 
                    climits[0] = np.min(var1) 
                    climits[1] = np.max(var1) 
                    #print('climits',climits)
                if j == 'u' or j == 'v' or j == 'w' or j == 'melt*_xy':
                    climits[0] = -1.*max(abs(climits[0]),abs(climits[1]))
                    climits[1] =     max(abs(climits[0]),abs(climits[1]))
            plt.clim(climits)
            
            #plt.title(runlabel + ' ' + t_label + ' = ' + str(int(tval)),
            #          fontsize = fs)

            plt.tight_layout()
            print(filenamesave)
            plt.savefig(outputdir+filenamesave)
            plt.close()
                    
    data1.close()
    if runcmp:
        datacmp.close()

#------------------------------------------------------------------------------
# PLOT_PR_TS
#   creates a TS plot from profile data colored by depth
# Inputs:
#   filedir  cell array of length 2 containing paths to 2 run directories
#   runname  cell array of length 2 containing shorthand labels for runs
#   teval    vector of times to plot field in hours (optional). Default is the 
#            first and last model output times. 
#   zlim     axis limits for depth (optional). Defaults to full model depth. 
#   zint     depth interval at which to plot velocity vectors (optional). 
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#   coupled  logical (optional). If true, plot both atmosphere and ocean 
#            velocity vectors
#------------------------------------------------------------------------------
def plot_TS(filedir,runname,teval = [-9999.], tunits = 'hr',
            zval = [9999.,9999.],zint=-10.,data_type='pr',
            outputdir = [], printformat = 'png',coupled = False):

    data1 = load_data(filedir,data_type=data_type)
    
    # set color = colorVal
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=zval[0], vmax=zval[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    for i in teval:
        fig = plt.figure()
        ax = fig.add_subplot(111)
   
        pt,y_axis_label = extract_var(data1,'pt',data_type=data_type,
                                         tval=[i,i],zval=zval)
        sa,x_axis_label = extract_var(data1,'sa',data_type=data_type,
                                         tval=[i,i],zval=zval)
            
        tval,t_label = extract_var(data1,'time',data_type='ts',tval=[i,i],keep='t',
                                   ops=tunits)
        z,cbar_label = extract_var(data1,'z',data_type=data_type,zval=zval)

        ax.scatter(sa,pt, s=ms, c=z,
                   cmap = cm, norm = cNorm, marker = 'o')
        ax.set_xlabel(x_axis_label,fontsize = fs)
        ax.set_ylabel(y_axis_label,fontsize = fs)
        
        cbar = fig.colorbar(scalarMap,ticks=np.arange(zval[0],zval[-1],zint*5))
        cbar.set_label(cbar_label)
        
        name = 'TS_z_scatter_'
        filenamesave = name + str(int(tval)) + tunits + '.' + printformat
        print(filedir + filenamesave)
        plt.savefig(filedir + filenamesave, format = printformat, bbox_inches="tight")
        plt.close()
    
    data1.close()
        
    return

#------------------------------------------------------------------------------
# PLOT_TSERIES
#   plot run variable vs. time
# Inputs:
#   filedir  path to run directory
#   runname  shorthand label for run
#   plotvar  cell array of variable names to be plotted in separate figures 
#   tlim     axis limits of time axis (optional). Vector of length 2.
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#------------------------------------------------------------------------------
def plot_tseries(filedir, runname, plotvar, 
                 tlim=[9999.,9999.],tav=0.,tunits='hr',tshade=[9999,9999],
                 ylim=[-9999,-9999], runlabel=[''], 
                 leglocation='best', legtitle='',legvar='', plot_legend=True,
                 col=col, linestyle=['-'], linewidth=pv.lw1, marker=mk, figsize=(6.4,4.8),
                 outputdir=[], overwrite=True, printformat='png'):
    
    if runlabel[0] == '':
        runlabel = runname
    print(runlabel)
    if len(linestyle) < len(filedir):
        linestyle = [linestyle[0] for i in filedir]

    for i in plotvar:
        
        if i in pv.varsname:
            varname = pv.varsvars[pv.varsname.index(i)]
            name = i + '_'
        else:
            varname = [i]
            name = pv.varname[pv.varlist.index(i)] + '_' 
        
        #name = pv.varname[pv.varlist.index(varname[0])] + '_'
        if len(runname) > 1:
            name = name + 'cmp_' + runname[1] + '_'
        if len(runname) > 2:
            name = name + runname[2] + '_'
        if tav > 0.:
            name = name + 'tav' + tav +'_'
        if tlim[0] != -9999.:
            name = name + 'tlim' + str(int(tlim[0])) + '-' + str(int(tlim[1])) +'_'
        filenamesave = name + 't'+'.' + printformat
        if len(outputdir) == 0:
            outputdir = filedir[0]
        if os.path.exists(outputdir + filenamesave) and not overwrite:
            print('skipping existing file: '+ outputdir + filenamesave)
            continue 
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)

        ymin = -9999
        ymax = -9999
        for idx,k in enumerate(filedir):
            data1 = load_data(k,data_type='ts')
            t,x_axis_label = extract_var(data1,'time',data_type = 'ts',
                                            tunits=tunits,tval=tlim,
                                            keep='t')
            if len(t) == 0:
                continue;
            
            for jidx,j in enumerate(varname):
               
               if len(filedir)>1:
                   if legvar != '':
                       ln_label = r'$'+varlabel[varlist.index(legvar)]+r'='+runlabel[idx]+r'$'
                   else:
                       ln_label = r'$' + runlabel[idx] + '\: $' 
               else:
                   ln_label = ''
               
               var,var_label = extract_var(data1,j,data_type = 'ts',
                                  tunits=tunits,tval=tlim,
                                  grid=pv.vartype[pv.varlist.index(j)],
                                  data_dir=k)
               if len(varname) > 1:
                   y_axis_label = r'$'+pv.vars_axis_label[pv.varsname.index(i)]+r'$'
                   linestyle[idx] = pv.vars_ls[pv.varsname.index(i)][jidx]
                   lw = pv.vars_lw[pv.varsname.index(i)][jidx]
                   ln_label = ln_label + var_label
               else:
                   y_axis_label = var_label
                   linestyle[idx] = '-'
                   lw = pv.lw1

               ax.plot(t,var,label = ln_label,
                       marker = marker,markersize = 0.5,
                       linestyle=linestyle[idx],linewidth=lw,
                       color = col[idx])

        if plot_legend:
            ax.legend(loc=leglocation, bbox_to_anchor = (1,1), title=legtitle)
        ax.set_xlabel(x_axis_label,fontsize = fs)
        if ylim[0] != -9999.:
           ax.set_ylim(ylim)
        if tlim[0] != 9999.:
           ax.set_xlim(tlim)
        if tshade[0] != 9999.:
           #ylim = ax.get_ylim()
           ax.fill([tshade[0],tshade[0],tshade[1],tshade[1]],
                   [ylim[0],ylim[1],ylim[1],ylim[0]],
                   'k',alpha = 0.2)#,lineStyle='None')
        ax.set_ylabel(y_axis_label,fontsize = fs)
        plt.yscale(pv.varscale[pv.varlist.index(varname[0])])
        
        print('save figure: '+outputdir + filenamesave)
        plt.savefig(outputdir + filenamesave, bbox_inches="tight")
        plt.close()
    
        data1.close()

#------------------------------------------------------------------------------
# PLOT_TSERIES_CROSS
#   creates a 
# Inputs:
#   filedir  cell array of length n containing paths to run directories
#   runname  cell array of length n containing shorthand labels for runs
#   teval    vector of times to plot field in hours (optional). Default is the 
#            first and last model output times. 
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#------------------------------------------------------------------------------
def plot_tseries_cross(filedir, runname, runlabel,
                       plotvar = ['us','melt'], data_type = ['ts','ts'],
                       teval = [-9999.,-9999.], tunits='hr', tav=0.,
                       zeval = [9999.,9999.],
                       figsize = (6.4,4.8), plot_cycles=False,
                       plot_jenkins=False, fit_power_law=False,
                       plot_legend = True, legtitle='',
                       x_ratio = False,
                       color_by_time = False, col = col, 
                       overwrite=False, outputdir = '', printformat = 'png'):
        
    name = pv.varname[pv.varlist.index(plotvar[1])] + '_' + pv.varname[pv.varlist.index(plotvar[0])]
    if len(runname) > 1:
        name = name + '_cmp_' + runname[1]
    if len(runname) > 2:
        name = name  + '_' + runname[2]
    if teval[0] != -9999.:
        name = name + '_tlim' + str(int(teval[0])) + '-' + str(int(teval[1]))
    if tav > 0.:
        name = name + '_tav' + str(int(tav))
    if plot_cycles:
        name = name + '_tcycles'
    if fit_power_law:
        name = name + '_fit'
    filenamesave = name + '.' + printformat
    if outputdir == '':
       outputdir = filedir[0]
    print(outputdir + filenamesave)
    if os.path.exists(outputdir + filenamesave) and not overwrite:
        print(outputdir + filenamesave, 'exists')
        return

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    data_all = np.zeros((2,len(filedir))) 
    for i,diri in enumerate(filedir):
        y_data = load_data(diri,data_type=data_type[1])
        if data_type[1] != 'parameter':
            slice_obj_input,varaxes = slice_var(y_data,'time',data_type=data_type[1],
                                                tunits=tunits,tval=teval)
        y_var,y_axis_label = extract_var(y_data,plotvar[1],
                                         slice_obj=slice_obj_input,
                                         data_type=data_type[1],data_dir=diri,
                                         tval=teval,tav=tav,zval=zeval)
        x_var_2,_ = extract_var(y_data,'u*',
                                slice_obj=slice_obj_input,
                                data_type=data_type[1],data_dir=diri,
                                tval=teval,tav=tav,zval=zeval)
        x_data = load_data(diri,data_type=data_type[0])
        if data_type[0] != 'parameter':
            slice_obj_x,varaxes = slice_var(x_data,'e*',data_type=data_type[0],
                                            tunits=tunits,tval=teval,zval=zeval)
            #slice_obj_x,varaxes = slice_var(x_data,'time',data_type=data_type[0],
            #                                tunits=tunits,tval=teval,zval=zeval)
        else:
            slice_obj_x = []
        x_var,x_axis_label = extract_var(x_data,plotvar[0],
                                         slice_obj=slice_obj_x,
                                         data_type=data_type[0],data_dir=diri,
                                         tval=teval,tav=tav,zval=zeval)
        if x_ratio:
            x_var = np.divide(x_var,x_var_2)
            x_axis_label = r'$\textrm{TKE}_{-2m}/u_*$'
        data_all[0,i] = x_var[0]
        data_all[1,i] = y_var[0]
        if color_by_time:
            t,cbar_label = extract_var(y_data,'time',data_type=data_type[1],tunits=tunits,tval=teval,
                                       keep='t',slice_obj=slice_obj_input)
            jet = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=np.min(t), vmax=np.max(t))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        
            ax.scatter(x_var, y_var, s=ms, c=t,
                       cmap=cm, norm=cNorm, marker='.')
            cbar = fig.colorbar(scalarMap)
            cbar.set_label(cbar_label)
        else:
            ax.plot(x_var, y_var, '.',
                    label = r'$'+runlabel[i]+r'$', c=col[i], 
                    markersize=10)
            if plot_legend:
                plt.legend(title=legtitle)
            if plotvar[0] == 'u*' and data_type[1] == 'pr':
                u,_ = extract_var(y_data,'U',data_type='pr',zval=[-2.,-2.],tval=[37.,50.])
                us = (0.003)**(0.5)*np.mean(u)
                ax.plot(us, y_var, 'o', label='', c=col[i], markersize=5, fillstyle='none')
        if plot_cycles:
            for j in np.arange(1,4):
                if data_type[0] != 'parameter':
                    x_var,_ = extract_var(y_data,plotvar[1],
                                          slice_obj=slice_obj_input,
                                          data_type=data_type[1],data_dir=diri,
                                          tval=teval-tav*j,tav=tav)
                add_var,_ = extract_var(y_data,plotvar[1],
                                        slice_obj=slice_obj_input,
                                        data_type=data_type[1],data_dir=diri,
                                        tval=teval-tav*j,tav=tav)
                ax.plot(x_var, add_var, '.', label='', c=col[i], 
                        #alpha = 0.8, 
                        markersize=10-2*j)

    if plot_jenkins:
        xmin,xmax = ax.get_xlim()
        ax.plot([xmin,xmax],[0.011,0.011],'--',c='green',linewidth=lw1)
        ax.set_xlim([xmin,xmax])
    if fit_power_law:
        print('power law fit is currently disabled')
        return
        pars, cov = curve_fit(f=power_law, xdata=data_all[0,:], ydata=data_all[1,:], 
                              p0=[0, 1], bounds=(-np.inf, np.inf))
        pars_set, cov_set = curve_fit(f=power_law_set, xdata=data_all[0,:], ydata=data_all[1,:], 
                                      p0=[0], bounds=(-np.inf, np.inf))
        # Get the standard deviations of the parameters (square roots of the # diagonal of the covariance)
        stdevs = np.sqrt(np.diag(cov))
        stdevs_set = np.sqrt(np.diag(cov))
        ax.plot(data_all[0,:],power_law(data_all[0,:],pars[0],pars[1]),'--k')
        # Calculate the residuals
        res = data_all[1,:] - power_law(data_all[0,:],pars[0],pars[1])
        res_set = data_all[1,:] - power_law_set(data_all[0,:],pars_set[0])
        r2_fit = rsquared(data_all[1,:],power_law(data_all[0,:],pars[0],pars[1]))
        r2_set = rsquared(data_all[1,:],power_law_set(data_all[0,:],pars_set[0]))
        print(r2_fit,pars)
        print(r2_set,pars_set)

    ax.set_xlabel(x_axis_label,fontsize = fs)
    ax.set_ylabel(y_axis_label,fontsize = fs)
    
    plt.savefig(outputdir + filenamesave, bbox_inches="tight")
    plt.close()
    
    x_data.close()
    y_data.close()

#------------------------------------------------------------------------------
# PLOT_UV_VECTOR
#   creates a quiver of velocity vector from profile data as a function of 
#   model depth
# Inputs:
#   filedir  cell array of length 2 containing paths to 2 run directories
#   runname  cell array of length 2 containing shorthand labels for runs
#   teval    vector of times to plot field in hours (optional). Default is the 
#            first and last model output times. 
#   zlim     axis limits for depth (optional). Defaults to full model depth. 
#   zint     depth interval at which to plot velocity vectors (optional). 
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#   coupled  logical (optional). If true, plot both atmosphere and ocean 
#            velocity vectors
#------------------------------------------------------------------------------
def plot_uv_vector(filedir,runname,runlabel = [''],
                   teval = 9999., tunits = 'hr',aspect='',
                   zval = [9999.,9999.],zint = 10.,zcolor = False,
                   legtitle = '', from3d=False,av=False,
                   outputdir = [], printformat = 'png',coupled = False):
    
    if runlabel[0] == '':
        runlabel = runname
    
    if from3d:
        data_type = '3d'
        av = True
    else:
        data_type = 'pr'
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(0, 0, marker = 'o', linestyle = 'none', color = 'none')
            
    for idx,i in enumerate(filedir):
        data1 = load_data(i,data_type=data_type,av=av)
        if coupled:
            data2 = load_data(filedir,coupled=True)
            var2,x_axis_label2 = extract_var(data2,j,tval=[teval,teval],zval=zval)
            z2,temp = extract_var(data1,'z',data_type=data_type,zval=zval)
        
        if len(filedir)>1:
            ln_label = r'$' + runlabel[idx] + '\: $' 
        else:
            ln_label = ''
            
        if teval == 9999.:
            teval = 1000.# this might not be necessary
        if zval[0] == -9999.:
            z,temp = extract_var(data1,'z',data_type=data_type)
            zval = [np.min(z),0]
        if zval[0] == zval[1]:
            zval[1] = zval[0]+zint
        
        tval,t_label = extract_var(data1,'time',tunits=tunits,
                                      data_type='ts',tval=[teval,teval])
        if idx == 0:
            filenamesave = 'u_vector_z_scatter' 
            if len(filedir) > 1:
                filenamesave = filenamesave + '_cmp_' + filedir[1]
            filenamesave = filenamesave + '_' + str(int(tval)) + tunits 
            if from3d:
                filenamesave = filenamesave + '_3dav'
            filenamesave = filenamesave + '.' + printformat # add zeval?
            outputdir = filedir[0] 
            if os.path.exists(outputdir + filenamesave) and not overwrite:
                print('skipping existing file: '+ outputdir + filenamesave)
                return

        if zcolor:
            jet = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=zval[0], vmax=zval[1])
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        else:
            #colorVal = np.zeros((4,))
            colorVal = col[idx]

        for j in np.arange(zval[0],zval[1],zint):
            zeval,cbar_label = extract_var(data1,'z',data_type=data_type,
                                              keep='z',zval=[j,j])
            u,x_axis_label   = extract_var(data1,'u',data_type=data_type,keep = 'z',
                                              zval=[j,j],tval=[teval,teval])
            v,y_axis_label   = extract_var(data1,'v',data_type=data_type,keep = 'z',
                                              zval=[j,j],tval=[teval,teval])
            if from3d:# consider replacing by ops=mean
                u = np.mean(u)
                v = np.mean(v)
        
            if zcolor:
               colorVal = scalarMap.to_rgba(j)
               ln_label = cbar_label + ' = ' + str(int(zeval))
          
            x = np.zeros((np.shape(u)))
            y = np.zeros((np.shape(u)))
            # plot points
            #ax.plot(u, v, label = ln_label,
            #        marker = '.', linestyle = 'none', color = colorVal)
            # plot vector
            #ax.plot([0,u],[0,v], label = ln_label,
            #        marker = mk,linestyle = ls, color = colorVal)
            ax.quiver(x,y,u,v,
                      angles='xy', scale_units='xy', scale=1,
                      color = colorVal, label = ln_label, linestyle = ls)
        
        if zcolor:
           cbar = fig.colorbar(scalarMap,ticks=np.arange(zval[0],zval[1],zint*2))
           cbar.set_label(cbar_label)
        elif len(filedir) > 1:
           ax.legend(loc=legloc, bbox_to_anchor=bboxanchor, title=legtitle)
        ax.set_xlabel(x_axis_label,fontsize = fs)
        ax.set_ylabel(y_axis_label,fontsize = fs)
        if aspect == 'equal':
            ax.set_aspect('equal','box')

    filenamesave = 'u_vector_z_scatter_' + str(int(tval)) + tunits 
    if from3d:
        filenamesave = filenamesave + '_3dav'
    filenamesave = filenamesave + '.' + printformat # add zeval?
    outputdir = filedir[0] 
    print(outputdir + filenamesave)
    plt.savefig(outputdir + filenamesave, bbox_inches="tight")
    plt.close()
    
    data1.close()

#------------------------------------------------------------------------------
# PLOT_PR
#   creates a plot of PALM variable field in YZ plane from 3D output data at user-
#   specified X
# Inputs:
#   filedir  path to run directory
#   runname  shorthand label for run
#   plotvar  cell array of variable names to be plotted 
#   teval    column 1 is minimum time, column 2 is maximum time
#            row corresponds to each run
#            vector of times to plot field in hours (optional). Default is the 
#            first and last model output times. 
#   tall     if true, plot all times between teval[i,0] and teval[i,1]
#   clim     axis limits for variable in plotvar (optional). 
#            Currently the same limits are applied to every variable in plotvar
#   zlim     axis limits for depth (optional). Defaults to full model depth. 
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#   coupled  logical (optional). if true, plot atmosphere and ocean profiles 
#            on the same axis
#   norm     string. if 'geostrophic', normalize by geostrophic velocity. 
#            if 'geostrophic_ice', normalize by geostrophic velocity for sloped
#            ice as in Jenkins 2016. 
#------------------------------------------------------------------------------
def plot_pr(filedir, runname, plotvar, 
            teval = [-9999.,-9999.], tall = False, tunits = 'hr', tav = 0.,
            xlim = [-9999.,-9999.], zlim = [9999.,9999.], 
            xscale_input = [], xscale_label = '', xscale = '', zscale = '', 
            data_type='pr',av=False,
            col=col, linestyle = ['-'], runlabel = [''], ops = [],
            plot_legend=True, legtitle = '',legvar = '', leglocation = legloc,
            show_boundary_value=True,hide_boundary_value_text = False,
            marker=mk, coupled = False,runcmp=False,figsize = (6.4,4.8),
            outputdir = [], printformat = 'png', overwrite=True, write_to_file = False 
            ):
    clim = xlim
    if xscale_input == []:
       xscale_input = np.ones((len(filedir)))

    if runlabel[0] == '':
        runlabel = runname
    if len(filedir) > 1:
        runcmp = True
    if len(linestyle) < len(filedir):
        linestyle = ['-' for i in filedir]
    if len(ops) < len(plotvar):
        ops = ['' for i in plotvar]
    if teval[0] == -9999.:
        if tall:
            # default values return full range of time should this be [9999.,9999.]? TODO
            tlim = [0.,9999.]
        else:
            #choose maximum time value
            tlim = [1000.,1000.]
    elif teval[0] != teval[1] and not tall:
        tlim = [teval[0],teval[0]]
    else:
        tlim = teval
    for i in plotvar:
        if i in pv.varsname:
            varname = pv.varsvars[pv.varsname.index(i)]
            name = i + ops[plotvar.index(i)] + '_' 
        else:
            varname = [i]
            name = pv.varname[pv.varlist.index(i)] + ops[plotvar.index(i)] + '_' 
        #print(i,varname)
        
        if len(runname) > 1:
            name = name + 'cmp_' + runname[1] + '_' 
        if len(runname) > 2:
            name = name + runname[2] + '_'
        if not tall:
            name = name + str(int(tlim[0]))+tunits+'_'
        else:
            name = name + 'tlim' + str(int(tlim[0])) + '-' + str(int(tlim[1])) +'_'
        if tav > 0:
            name = name + 'tav' + str(tav)+ '_'
        if data_type == '3d':
            name = name + '3dav_'
        if zlim[0] != 9999.:
            name = name + str(int(abs(zlim[0]))) + 'zmax_'
        if xscale_label != '':
            name = name + 'scale_'
        if len(outputdir) == 0:
            outputdir = filedir[0]
        filenamesave = name + 'z_profile.' + printformat
        if write_to_file:
            df = open(name + 'z_profile.txt','a+')
            col_headings=['z']
        if os.path.exists(outputdir + filenamesave) and not overwrite:
            print('skipping existing file: '+ outputdir + filenamesave)
            continue 
        print('generating file: '+ outputdir + filenamesave)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        cmin = 9999
        cmax = -9999
       
        for idx,k in enumerate(filedir):
            data1 = load_data(k,data_type=data_type,av=av)
            if coupled:
                data2 = load_data(filedir,coupled=True)
                var2,x_axis_label2 = extract_var(data2,j);
            
            #slice_obj_input,varaxes = slice_var(data1,varname[0],data_type='pr',
            #                              tunits=tunits, tval=tlim, zval=zlim)
            #print(f'plot_pr: slice_obj={slice_obj_input}')
            #print(f'plot_pr: varaxes={varaxes}')
            t,t_label= extract_var(data1,'time',data_type='ts',
                                   #slice_obj=[slice_obj_input[varaxes.index('t')]],
                                   keep='t',tunits=tunits,tval=tlim)
            if tall:
                jet = cm = plt.get_cmap('jet') 
                cNorm  = colors.Normalize(vmin=min(t), vmax=max(t))
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

            # load z 
            z,y_axis_label = extract_var(data1,'z',zval=zlim,
                             #slice_obj=[slice_obj_input[varaxes.index('z')]],
                             grid=pv.vartype[pv.varlist.index(varname[0])])
            if zscale == 'Ekman':
               zE,_ = extract_var(data1,'zE',zval=zlim)
               #zE = derived_var(data1,'zE',slice_obj_input)
               z = np.divide(z,zE)
               #print('z=',z)
               y_axis_label = r'$z/d_E$'
            
            if coupled:
                z2,y_label = extract_var(data2,'z',zval=zlim,
                             grid=pv.vartype[pv.varlist.index(varname[0])])
            
            for jidx,j in enumerate(varname):

                if len(filedir)>1:
                    if legvar != '':
                        ln_label = r'$'+varlabel[varlist.index(legvar)]+r'='+runlabel[idx]+r'$'
                    else:
                        ln_label = r'$' + runlabel[idx] + '\: $' 
                else:
                    ln_label = ''
                
                # load z only if different variables are on differnt grids
                var1,var_label = extract_var(data1, j, zval=zlim,
                                             tunits=tunits, tval=tlim, keep='t', tav = tav,
                                             data_type=data_type,ops=ops[plotvar.index(i)],
                                             #slice_obj=slice_obj_input,
                                             data_dir=k)
                if pv.vartype[pv.varlist.index(j)] != pv.vartype[pv.varlist.index(varname[0])]: 
                    z,y_axis_label = extract_var(data1,'z',zval=zlim,
                                                 grid=pv.vartype[pv.varlist.index(j)])
                    if coupled:
                        z2,y_label = extract_var(data2,'z',zval=zlim,
                                                 grid=pv.vartype[pv.varlist.index(j)])
                
                if coupled:
                    var1, = extract_var(data2,j,tval=tlim,keep='t',
                                           data_type=data_type,ops=ops[plotvar.index(i)],
                                           data_dir=k)
                if len(varname) > 1:
                    ls = pv.vars_ls[pv.varsname.index(i)][jidx]
                    lw = pv.vars_lw[pv.varsname.index(i)][jidx]
                    ln_label = ln_label + var_label
                else:
                    ls = linestyle[idx]
                    lw = lw1
               
                if data_type=='3d':
                    var1 = np.mean(np.mean(var1,axis=1),axis(1))
                
                if np.shape(var1)[1] == 0:
                    print(k+'skip because '+j+' is empty')
                    continue

                for tidx,tval in enumerate(t):
                    if tall:
                        colorVal = scalarMap.to_rgba(tval)
                    else:
                        colorVal = col[idx]
                    #print('var1=',np.divide(var1[tidx,:],xscale_input[idx]))
                    #print('shape(var1)=',np.shape(var1))
                    #print('np.mean(var1)=',np.mean(var1))
                    ln, = ax.plot(np.divide(var1[tidx,:],xscale_input[idx]), z, 
                                  label = ln_label, linewidth=lw,
                                  marker = marker,linestyle = ls,
                                  color=colorVal)
            
                if coupled:
                    ax.plot(var2[tidx,:], z2, linewidth=lw,
                            marker = marker,linestyle = ls,color=colorVal)
                if not show_boundary_value:
                    cmin = min(cmin,np.min(var1[tidx,:-1]))
                    cmax = max(cmax,np.max(var1[tidx,:-1]))
                    if not hide_boundary_value_text:
                        ax.text(clim[0] + 0.05*(clim[1] - clim[0]),
                              z[-2] + 0.5*(  z[-1] - z[-2]  ),
                            '{:3.2f}'.format(var1[tidx,-1]),
                            fontsize = fs1,color=colorVal)
                #if write_to_file:
                #   for k in filedir:
                #      col_headings.append(k+varname)
                #df.writerow()
            data1.close()
            
        if not tall:
           if plot_legend:
              #ax.legend(loc=legloc, bbox_to_anchor=bboxanchor, title=legtitle)
              #ax.legend(loc='lower right', title=legtitle)
              ax.legend(title=legtitle)
        if len(varname) > 1:
            var_label = r'$'+pv.vars_axis_label[pv.varsname.index(i)]+r'$'
        if xscale_label != '':
            var_label1,var_label2 = var_label.split('(')
            var_label = var_label1 + xscale_label
            #var_label = var_label1 + xscale_label + var_label2

        if i == 'velocity':
            clim = ([-1.*max(abs(cmin),abs(cmax)),
                         max(abs(cmin),abs(cmax))])
        #else:
        #    clim = [cmin,cmax]
        if not show_boundary_value:
            clim[1] += 0.05*(clim[1] - clim[0])
        if tall:
            cbar = fig.colorbar(scalarMap)#,ticks=t)
            cbar.set_label(t_label,fontsize=fs)
        if xscale_label != '':
           ax.set_ylim([-25,0])
        if abs(zlim[0]) != 9999.:
            ax.set_ylim(zlim)
        if abs(clim[0]) != 9999.:
            ax.set_xlim(clim)

        ax.set_xlabel(var_label)
        ax.set_ylabel(y_axis_label)
        if xscale != '':
            plt.xscale(xscale)
        else:
            plt.xscale(pv.varscale[pv.varlist.index(varname[0])])
            #plt.xscale(pv.varscale[pv.varlist.index(i)])
        if not runcmp and not tall:
            plt.title(runlabel[0] + ' ' + t_label + ' = ' + str(int(t)))
        
        print('saving file: '+ outputdir + filenamesave)
        plt.savefig(outputdir + filenamesave, bbox_inches="tight")
        plt.close()
        clim = xlim
    return

# either designate multiple times or multiple directories
def plot_hovmoller(filedir, runname, plotvar,tunits='hr', 
                   runlabel = [''], 
                   contour_var = '', contour_val = 0, contour_col = 'thistle',
                   plot_legend=True, plot_BL=False,
                   zlim = [-9999,-9999], clim = [-9999.,-9999.],tlim = [-9999.,-9999.],
                   figsize = (6.4,4.8), outputdir = [], printformat = 'png', overwrite=False):
    if runlabel[0] == '':
        runlabel = runname
    for j in plotvar:

        name = pv.varname[pv.varlist.index(j)]# + ops[plotvar.index(j)] 
        #if len(runname) > 1:
        #    name = name + '_cmp_' + runname[1]
        if len(runname) > 2:
            name = name + '_' + runname[2]
        if tlim[0] != -9999.:
            name = name + '_tlim' + str(int(tlim[0]))+'-'+ str(int(tlim[1]))
        else:
            name = name + '_t'
        if zlim[0] != -9999.:
            name = name + '_' + str(int(abs(zlim[0]))) + 'zmax'
            #name = name + '_z{0.03d}-{1.03d}'.format(zlim[0],zlim[1]) + 'm'
        if plot_BL:
            name = name + '_BL'
        if contour_var != '':
            name = name + '_contour_' + contour_var
        #print(name)
        filenamesave = name + '_hovmoller.' + printformat
        
        for idx,k in enumerate(filedir):
            outputdir = k
            if os.path.exists(outputdir + filenamesave) and not overwrite:
                print('skipping existing file: '+ outputdir + filenamesave)
                continue

            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(111)
        
            ln_label = runlabel[idx]
            data1 = load_data(k)
            var1,c_axis_label = extract_var(data1,j,#ops=ops[plotvar.index(j)],
                                            data_type = 'pr',tunits=tunits, 
                                            tval=tlim, zval=zlim, data_dir=k)
            zu,y_axis_label = extract_var(data1,'z',data_type = 'z',
                                          grid='u', zval=zlim)
            zw,_            = extract_var(data1,'z',data_type = 'z',
                                          grid='w', zval=zlim)
            t,x_axis_label  = extract_var(data1,'time', data_type = 'ts',
                                          tunits=tunits, tval=tlim)
            if (pv.varscale[pv.varlist.index(j)]== 'log'):
                var1 = np.log10(var1)
                c_axis_label = r'$\log_{10}$'+c_axis_label

            if contour_var != '':
                if contour_var != j:
                    var2,_ = extract_var(data1,contour_var,#ops=ops[plotvar.index(j)],
                                         data_type = 'pr',tunits=tunits, 
                                         tval=tlim, zval=zlim, data_dir=k)
                else:
                    var2 = var1
                contour_z = np.ones((len(t)))*np.min(zw)
                for idx_t,ti in enumerate(t):
                    idx_z = 0
                    while var2[idx_t,idx_z] < contour_val and idx_z < len(zw)-1:
                        contour_z[idx_t] = zw[idx_z]
                        idx_z += 1
                ax.plot(np.add(t,0.5),contour_z,'--',color=contour_col,linewidth=lw1)
            if plot_BL:
                BL_z = np.zeros((len(t)))
                pt,_ = extract_var(data1,'pt',#ops=ops[plotvar.index(j)],
                                   data_type = 'pr',tunits=tunits, 
                                   tval=tlim, zval=zlim, data_dir=k)
                dT,_ = extract_var(data1,'thermal_driving_infty',
                                   data_type = 'parameter',tunits=tunits, 
                                   tval=tlim, zval=zlim, data_dir=k)
                pt_f = gsw.pt0_from_t(35,gsw.t_freezing(35,800,0),800)
                n_crit = 0.99
                pt_BL = n_crit * dT[0] + pt_f
                for idx_t,ti in enumerate(t):
                   BL_z[idx_t] = zw[np.argmin(np.abs(np.subtract(pt[idx_t,:],pt_BL)))]
                #print(pt_BL,dT+pt_f)
                #print(BL_z)
                ax.plot(np.add(t,0.5),BL_z,'-',color='k',linewidth=lw1)
            Z = np.shape((len(zu)+1,len(t)+1))# row is z,col is t 
            T = np.shape((len(zu)+1,len(t)+1))# row is z,col is t 
            zw = ma.append(zw[0]+(zw[0]-zw[1]),zw)
            t = ma.append(t,t[-1] + (t[-1]-t[-2]))
            Z,T = np.meshgrid(zw,t)
            
            if clim[0] == -9999:
                cNorm  = colors.Normalize(vmin=np.percentile(var1,5), vmax=np.percentile(var1,95))
            else:
                cNorm = colors.Normalize(vmin=clim[0], vmax=clim[1])
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pv.varcmap[pv.varlist.index(j)])
            pmesh = ax.pcolormesh(T,Z,var1,cmap = pv.varcmap[pv.varlist.index(j)], norm=cNorm)
            data1.close()
           
            if zlim[0] != -9999.:
                ax.set_ylim(zlim)
            if tlim[0] != -9999.:
                ax.set_xlim(tlim)
            
            #if len(filedir) > 1:
            #    ax.legend(loc = 9,bbox_to_anchor=(0.5, -0.15))
            if plot_legend:
                cbar = fig.colorbar(scalarMap)
                cbar.set_label(c_axis_label)
            ax.set_xlabel(x_axis_label,fontsize = fs)
            ax.set_ylabel(y_axis_label,fontsize = fs)
            
            print(outputdir + filenamesave)
            plt.savefig(outputdir + filenamesave, bbox_inches="tight")
            plt.close()

#------------------------------------------------------------------------------
# PLOT_TSERIES_FROM_PR
#   plot variable vs. time at a user-specified z-level or if 
#   formerly PLOT_TSERIES_ZLEVEL
# Inputs:
#   filedir  path to run directory
#   runname  shorthand label for run
#   plotvar  cell array of variable names to be plotted in separate figures 
#   zeval    depth at which to plot variable (optional).
             # TODO for w*2_max, assign to 9999.
#   clim     axis limits of variable axis (optional). Vector of length 2.
#   tlim     axis limits of time axis (optional). Vector of length 2.
#   output_dir filepath to save plot (optional). Defaults to current directory.
#   printformat plot file extension (optional). Defaults to png.
#------------------------------------------------------------------------------
def plot_tseries_from_pr(filedir, runname, plotvar,ops=[''], 
                         runlabel = [''], legtitle = '',legvar = '',
                         zeval = 0., clim = [-9999.,-9999.],
                         tlim = [9999.,9999.],tunits='hr',
                         col=col, linestyle = ['-'], outputdir = [], 
                         printformat = 'png', overwrite=False):
    if runlabel[0] == '':
        runlabel = runname
    if len(linestyle) < len(filedir):
        linestyle = ['-' for i in filedir]
    
    for j in plotvar:

        if len(outputdir) == 0:
            outputdir = filedir[0]
        if j in pv.varsname:
            varname = pv.varsvars[pv.varsname.index(j)]
            name = j + ops[plotvar.index(j)] + '_' 
        else:
            varname = [j]
            name = pv.varname[pv.varlist.index(j)] + ops[plotvar.index(j)] + '_' 
        if len(runname) > 1:
            name = name + '_cmp_' + runname[1]
        if len(runname) > 2:
            name = name + '_' + runname[2]
        if tlim[0] != 9999.:
            name = name + '_tlim' + str(int(tlim[0]))+'-'+ str(int(tlim[1]))
        else:
            name = name + '_t'
        if zeval != 9999.:
            name = name + '_z' + str(int(abs(zeval)))+'m'
            #name = name + 't_'+str(int(abs(z)))+'m'
        filenamesave = name + '.' + printformat
        if os.path.exists(outputdir + filenamesave) and not overwrite:
            print('skipping existing file: '+ outputdir + filenamesave)
            continue 
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for idx,k in enumerate(filedir):
            data1 = load_data(k)
            t,x_axis_label = extract_var(data1,'time', data_type = 'ts',
                                            tunits=tunits, tval=tlim)
            for iidx,i in enumerate(varname):
                var1,y_axis_label = extract_var(data1,i,ops=ops[plotvar.index(j)],
                                                   data_type = 'pr',tunits=tunits, 
                                                   tval=tlim, zval=[zeval,zeval])
                #print(np.min(var1),np.max(var1))
                if legvar != '':
                    ln_label = r'$'+varlabel[varlist.index(legvar)]+r'='+runlabel[idx]+r'$'
                else:
                    ln_label = runlabel[idx]
                if len(varname) > 1:
                    ls = pv.vars_ls[pv.varsname.index(j)][iidx]
                    lw = pv.vars_lw[pv.varsname.index(j)][iidx]
                    ln_label = ln_label + ', ' + y_axis_label
                else:
                    ls = linestyle[idx]
                    lw = lw1

                #if zeval != 9999.:
                #   ln_label += r' z='+str(int(zeval))+' m'
                
                ax.plot(t,var1,
                     label = ln_label,
                     marker = mk, markersize = ms, 
                     linestyle = ls, linewidth = lw,color=col[idx])
                
            data1.close()
        if len(varname) > 1:
            y_axis_label = r'$'+pv.vars_axis_label[pv.varsname.index(j)]+r'$'
           
        if clim[0] != -9999.:
            ax.set_ylim(clim)
        if tlim[0] != -9999.:
           ax.set_xlim(tlim)
        plt.yscale(pv.varscale[pv.varlist.index(j)])
        
        if len(filedir) > 1:
            ax.legend(loc = 9,bbox_to_anchor=(0.5, -0.15),title=legtitle)
        ax.set_xlabel(x_axis_label,fontsize = fs)
        ax.set_ylabel(y_axis_label,fontsize = fs)
        
        print(outputdir + filenamesave)
        plt.savefig(outputdir + filenamesave, bbox_inches="tight")
        plt.close()

def rsquared(y,yfit):
    s_tot = np.sum(np.square(y-np.mean(y)))
    s_res = np.sum(np.square(y-yfit))
    return 1-(s_res/s_tot)

def power_law(x,a,b):
    return a*np.power(x,b)

def power_law_set(x,a):
    b = 4/3
    #b = 1
    return a*np.power(x,b)
