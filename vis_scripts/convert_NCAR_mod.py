"""
Created on Wed Jan 30 13:09:33 2019

@author: cbegeman
"""
from netCDF4 import Dataset
import numpy as np

def convert_NCAR(LES1):
    
    s_hr = 3600
    LES1.variables['time'] = np.multiply(LES1.variables['Time'],s_hr)
    LES1.variables['zu'] = LES1.variables['zt']
    LES1.variables['zw'] = LES1.variables['zm']
    LES1.variables['pt'] = np.add(LES1.variables['mean_temperature'],273.15)
    LES1.variables['sa'] = LES1.variables['mean_salinity']
    LES1.variables['u'] = LES1.variables['mean_U']
    LES1.variables['v'] = LES1.variables['mean_V']
    LES1.variables['w'] = LES1.variables['mean_W']
    LES1.variables['u*2'] = LES1.variables['lesResolved_uu']
    LES1.variables['v*2'] = LES1.variables['lesResolved_vv']
    LES1.variables['w*2'] = LES1.variables['lesResolved_ww']
    LES1.variables['w*u*'] = LES1.variables['lesResolveduwFlux'] #1e-5
    LES1.variables['w"u"'] = LES1.variables['subGridScheme_uw'] #1e-8
    LES1.variables['w*v*'] = LES1.variables['lesResolvedvwFlux']
    LES1.variables['w"v"'] = LES1.variables['subGridScheme_vw']
    LES1.variables['w*pt*'] = LES1.variables['lesResolvedwtFlux'] #1e-5
    LES1.variables['w"pt"'] = LES1.variables['subGridScheme_wt']
    LES1.variables['w*sa*'] = LES1.variables['lesResolvedwsFlux']
    LES1.variables['w"sa"'] = LES1.variables['subGridScheme_ws']
    LES1.variables['e'] = LES1.variables['sgsE']
    LES1.variables['e*'] = 0.5*(LES1.variables['u*2'][:] + 
                                  LES1.variables['v*2'][:] +
                                  LES1.variables['w*2'][:] )

    return LES1
