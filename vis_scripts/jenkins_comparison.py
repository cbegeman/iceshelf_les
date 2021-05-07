#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 07:52:34 2019

@author: cbegeman
"""
import math
import numpy as np

pt_far = -2.385
pt0 = [-2.505,-2.505,-2.5025,-2.5]
dT = [pt_far - i for i in pt0]
alpha = [1.e-2,1.e-1,5.e-1,1.]
ug_i = 17.5* np.multiply(np.sin(alpha),dT)
print(dT)
print(ug_i)
