# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from extract_var_palm import load_data,extract_var
import sys
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize2_wide,figsize2_box
#from cmp_thermal_driving_slope1 import *
#from cmp_slope_hidT import *

#sys.path.append('/Users/cbegeman/Software_files/my_python_code/')
teval = [10]
zmax = -64
filedir = '/lustre/scratch5/cbegeman/palm/jobs/test-chicoma-partialice-8/RUN_ifort.chicoma_hdf5_srun_test_partialice/'
palm.plot_3d_slice(filedir,'',['pt'],
                   teval = teval,
                   slice_dim = 'y',#clim=[-2.44,-2.39],
                   zval=[zmax,0],yval=[64,64],
                   figsize=figsize2_wide)
