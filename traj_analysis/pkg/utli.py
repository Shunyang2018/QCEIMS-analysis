#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 20:33:26 2019

@author: shunyang
"""
import re
import numpy as np



def process(coord):
    
    cartesian = np.zeros((len(coord)-1,3))              
    initial_line = re.split(r' +',coord[0]) # initial line contains several important info about the point
    atom_list = []
    
    for i in range(1,len(coord)):# extract out the coordinate of each atom and put it into a list
        
        line = re.split(r'[\s]',coord[i])
        atom_list.append(line[0])
        line = [x for x in line if x!= '']
        cartesian[i-1,0] = float(line[1])
        cartesian[i-1,1] = float(line[2])
        cartesian[i-1,2] = float(line[3])          
            
    return cartesian,initial_line,atom_list