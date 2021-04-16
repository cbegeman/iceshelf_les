"""
Created on Tue Jan 29 11:54:36 2019

@author: cbegeman
"""

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) or 1 # use start += 1 to find overlapping matches

def rm_palm_char(j):
    name = j + '_'
    
    ridx = list(find_all(j,'*'))
    if len(ridx) >= 1:
        name = j[0:ridx[0]] + '_res_'
    if len(ridx) == 2:
        name = name + j[ridx[0]+1:ridx[1]] + '_res_'
                                        
    sidx = list(find_all(j,'"'))
    if len(sidx) >= 1:
        name = j[0:sidx[0]] + '_sgs_'
    if len(sidx) == 2:
        name = name + j[sidx[0]+1:sidx[1]] + '_sgs_'
    
    return name;
