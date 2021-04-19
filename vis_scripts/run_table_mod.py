
import numpy as np
import os
import pysed as sed
import shutil
from namelist_default import ininame,inivalue
import gsw
import filecmp
import csv
import sys
import math

def load_rignot():
    filename = 'rignot_2002_dT_melt.txt'
    table_file = open(filename,'r')
    rd = csv.reader(table_file,delimiter=' ')
    glacier_year = table_file.readline().split(' ')
    melt = table_file.readline().split(' ')
    dT = table_file.readline().split(' ')
    
    print(np.shape(glacier_year))
    print(np.shape(melt))
    print(np.shape(dT))
    print(melt[-1])
    
    melt_mean = np.ones((len(melt)-1))
    melt_sd = np.ones((len(melt)-1))
    for idx,i in enumerate(melt_mean):
        dT[idx+1] = float(dT[idx+1])
        melt_mean[idx] = float(melt[idx+1][0:melt[idx+1].index('_')])
        melt_sd[idx] = float(melt[idx+1][melt[idx+1].index('_')+1:])
    return melt_mean,melt_sd,dT[1:]
def valfromline(line):
    varvalue = [0.,0.]
    string,value = line.split('=')
    varlist = value.split(',')
    for i,var in enumerate(varlist):
       if ( var.strip() == '.T.'):
           varvalue[i] = 1.
       elif ( var.strip() == '.F.'):
           varvalue[i] = 0.
       elif ( var.startswith("'") ):
           print('Other string\n')
           break
       elif var == '':
           break
       else:
          varvalue[i] = float(var.strip())
          break
    return [varvalue[i]]
def load_vals(filename):
   table_file = open(filename,'r')
   rd = csv.reader(table_file,delimiter=',')
   varname = table_file.readline().split(',')
   runname = []
   data = list(rd)
   numruns = len(data)
   derivedvar = ['Lx','Ly','Lz']
   numvars = len(data[1]) + len(derivedvar)
   runtable = np.zeros((numruns,numvars))
   varname[-1] = varname[-1][:-1]
   for var in derivedvar:
      varname.append(var)
   for i,row in enumerate(data):
      runname.append(row[varname.index('rundir')])
      for j,entry in enumerate(row):
         if j == 0 or not entry:
            row[j] = math.nan
         if j > 0 and type(row[j]) != float:
            row[j] = float(row[j])
      for var in derivedvar:
         if var[0] == 'L':
            L = row[varname.index('n'+var[1])] * row[varname.index('d'+var[1])]
            row.append(L)
      runtable[i,:] = row
   return varname,runname,runtable
def value_from_namelist(filename,var):
   datafile = open(filename+'/test_oceanml_p3d','r')
   for line in datafile:
       if var in line:
           return valfromline(line)
   return math.nan
