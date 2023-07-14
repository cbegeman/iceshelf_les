# -*- coding: utf-8 -*-
"""
Timeseries figure

"""

import sys
sys.path.append('../.')
from extract_var_palm import load_data,extract_var
import plot_palm_mod as palm
import numpy as np
from plot_param_palm import figsize3


filedir = '/lustre/scratch5/cbegeman/palm/jobs/test-chicoma-partialice-8/RUN_ifort.chicoma_hdf5_srun_test_partialice/'
palm.plot_tseries([filedir], [''], ['E*','u*'],
                  figsize=figsize3, plot_legend=False)
