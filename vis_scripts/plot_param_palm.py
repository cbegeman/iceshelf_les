#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 07:52:34 2019

@author: cbegeman
"""
# define plotting parameters
lw1 = 1
lw2 = 2
lw3 = 3 


ls = '-' #linestyle
#ls = 'None' #linestyle
mk = 'None' #marker
#mk = '.' #marker
ms = 0.25 #marker
fs = 16
fs1 = 14
legloc = 'upper left'
bboxanchor = (0.0,-0.25)
s_hr = 3600
logmin = 1.0e-20

col2_width = 12 #cm
col1_width = 8.3 #cm
figsize2_wide = (col2_width/2,2.5)
figsize2_box = (col2_width/2,3)
figsize1 = (col1_width,4)
figsize2 = (col2_width/2,col2_width/2)
figsize3 = (col2_width/3,col2_width/3)

# define default colors
col = ['#186086','#F5793A','#A95AA1','#85C0F9','#D14A42','#2ca25f','#99d8c9','#5F86DD']
#col = np.matrix('1 0 0;0 1 0;0 0 1')

