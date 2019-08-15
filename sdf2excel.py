#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:20:12 2019

@author: shunyang
transform sdf, ie mainlib.sdf to csv file! Watch out they are very big!
"""
import pandas as pd
import numpy as np


sdf = '/Users/shunyang/Box/qceims_input/code/test/test.sdf'


def read_blocks(file):
    block = []
    for line in file:
        if line.startswith('$$$$') and len(block)>0:
            yield block
            block = []
        block.append(line)
    yield block

def mass2list(mass):
    x=np.zeros([500,])# mass from 0 to 500
    
    for line in mass:
        temp = line.strip().split()

        try:    
#            print(int(temp[0]), float(temp[1]))
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
timer = 1
with open(sdf) as f:
    for block in read_blocks(f):
#        print(timer,'compounds')
        timer += 1
        n = len(block)
        flag = False
        for i in range(0,n):
            if '<NAME>' in block[i]:      
                NAME.append(block[i+1])
            if '<INCHIKEY>' in block[i]:
                INCHIKEY.append(block[i+1])
                flag = True
            if '<FORMULA>' in block[i]:
                FORMULA.append(block[i+1])
            if '<MASS SPECTRAL PEAKS>' in block[i]:
                mass = block[i+1:-1]
                test = mass2list(mass)
                MASS.append(test)
        if flag == False:
            print(timer, 'No-Inchikey')
            INCHIKEY.append('No-Inchikey')
INCHIKYEY = INCHIKEY.pop()           
print('NAME:',len(NAME),'INCHIKEY',len(INCHIKEY), 'FORMULA', len(FORMULA), 'MASS', len(MASS))               
                
df = pd.DataFrame({'NAME':NAME, 'INCHIKEY':INCHIKEY, 'FORMULA':FORMULA, 'MASS':MASS  })
#df['SHORT'] = df['INCHIKEY'].str.split('-', n=1,expand=True)
df.to_csv('test.csv')