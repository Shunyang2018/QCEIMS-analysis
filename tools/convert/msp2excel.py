#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:20:12 2019

@author: shunyang

We export MS from NIST by MSP file. To calculate simliarity, we want to use uniform format. So, we transform it
to csv file. To save the mass data, we introduce a 500 length list, whose index represents the m/z and value means the intensity.
 
"""
import pandas as pd
import numpy as np

file = '/Users/shunyang/Box/qceims_input/conformer-2-nonene/2-nonene.MSP'



def mass2list(mass):
    x=np.zeros([500,])# mass from 0 to 500
    
    for line in mass:
        temp = line.strip().split()

        try:    
            x[int(temp[0])] = float(temp[1])
        except IndexError:
            break
        except ValueError:
            break
    return x
        
        

NAME = []
INCHIKEY = []
FORMULA = []
MASS = []
timer = 0 # count if we get different number
with open(file) as f:
    block = f.readlines()
    block = [x.strip() for x in block] 
    timer += 1
    n = len(block)
    flag = False
    for i in range(0,n):
        if 'Name' in block[i]:      
            NAME.append(block[i].split(':'[1]))
        if 'InChIKey' in block[i]:
            INCHIKEY.append(block[i].split(':'[1]))
            flag = True
        if 'Formula' in block[i]:
            FORMULA.append(block[i].split(':'[1]))
        if 'Num Peaks:' in block[i]:
            mass = block[i+1:]
            test = mass2list(mass)
            MASS.append(test)
    if flag == False:
        print(timer, 'No-Inchikey')
        INCHIKEY.append('No-Inchikey')
INCHIKYEY = INCHIKEY.pop()           
print('NAME:',len(NAME),'INCHIKEY',len(INCHIKEY), 'FORMULA', len(FORMULA), 'MASS', len(MASS))               
# save to dataframe and .csv              
df = pd.DataFrame({'NAME':NAME, 'INCHIKEY':INCHIKEY, 'FORMULA':FORMULA, 'MASS':MASS  })
df.to_csv('2-nonene.csv')