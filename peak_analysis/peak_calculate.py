#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 12:27:53 2021

@author: shunyang
"""
import pickle
import pandas as pd
import numpy as np
import Daphnis.methods
from jdxtool import cos_sim, dot_sim, readjdx, peak, match,  trans




with open('index_insilico_ref.pkl','rb') as r:
    
    
    # read list of ins, ref spectra as [massindex, insilico, ref]
    spclist = pickle.load(r)
    
    
CM= np.empty(shape=[0,4]) # confusion matrix 
#examples:

ins = spclist[0][1]
ref = spclist[0][2]
dot = dot_sim(ins,ref) # modified score
cm = peak(ins,ref)  
CM = np.append(CM, cm, axis=0)