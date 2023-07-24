# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
sys.path.append('../.')
from extract_var_palm import load_data,extract_var
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize2

filedir = '/lustre/scratch5/cbegeman/palm/jobs/test-chicoma-partialice-8/RUN_ifort.chicoma_hdf5_srun_test_partialice/'
tprofile = 1
zmax = -9999 # For final figures, this should be replaced by 2/3 the domain size to exclude the sponge layer

palm.plot_pr([filedir], [''], 
             ['heatflux_z', 'saltflux_z',
              'momflux_u','momflux_v','momflux_z'],
             teval = [tprofile,tprofile], 
             zlim=[zmax,0],
             figsize = figsize2)
