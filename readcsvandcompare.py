#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 17:59:22 2019

@author: shunyang
"""


import numpy as np
import pandas as pd
import re

def valid_component_str(comp_str):
'''
When you use some other guys' codes, must make sure you understand it and double check it can give you the result consistent to your own program.
For this, like O and N, like the space in front of O and the lack of '0'.
Otherwise, you will take a long time to debug!!!


when you don't know where is the bug, go to the root!
'''
    cho_set = {'C', 'H','O', 'N'}
    residual1 = set(comp_str) - cho_set
    numerical_set = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'}
    residual2 = residual1 - numerical_set
    if len(residual2) == 0:
        return True
    else:
        return False

df = pd.read_excel('/Users/shunyang/Box/qceims_input/qceims_input.xlsx', header = 0, sheet_name = 'Sheet1')
df_des = pd.read_excel('/Users/shunyang/Box/qceims_input/NIST17-descriptors.xlsx',header=0, sheet_name=1)


#boxplot = df.boxplot(column=['Match'], by=['unsaturation '])
#df['inchikey'] == df_des['InchiKey']


df_des['Formula'] = df_des['Formula'].fillna('FFF')
test = df_des[df_des['Formula'].str.contains(r'C\dH\d')]

test = test[test['Formula'].apply(valid_component_str)]
print(test)
test.to_excel('clean-full.xlsx', sheet_name = 'test')
df = df.rename(columns={'short inchikey':'InchiKey-Short'})
df = df[0:235][:]
temp = test.copy()

temp = pd.merge(temp,df,  on='InchiKey-Short', how='inner')
#temp.to_excel('resultmerge.xlsx', sheet_name = 'test')

print(temp)
#T = []
#for i in df['InchiKey-Short']:
#    print(i)
#    t = temp[temp['InchiKey-Short'].str.contains(i)]
#    T.append(t)
#print(T)
    
       