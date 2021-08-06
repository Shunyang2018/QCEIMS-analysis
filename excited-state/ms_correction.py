#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 12:03:57 2021

@author: shunyang
"""
from jdxtool import cos_sim, dot_sim, readjdx
import pandas as pd
import numpy as np


def readblocks(file):
    block = []
    flag = 0
    for line in file:
#        if line.startswith('\n'):
#            print('not line')
#            flag = 0
#
#            block = []
        if flag == 1:
            block.append(line)
        if line.startswith('Num Peaks:'):
            flag = 1
    yield block



def readmsp(msp):
    with open(msp) as f:
    #another way to read line by line, which is useful in small flies
    #    test = f.readlines()
    #    print(test[0])
        for block in readblocks(f):
            x=np.zeros([699,])
            for line in block:
                line = line.strip('\n')
                line = line.split(';')
                for pair in line:
                    pair = pair.split()
                    if len(pair) == 2:
                        x[int(pair[0])] = int(pair[1])
    return x

#%%
path = '/Users/shunyang/project/xcited-state/qceims/MethylHistidine/'
msp = '/Users/shunyang/project/xcited-state/qceims/MethylHis.msp'#the experimental data
x = readmsp(msp)
file = '/Users/shunyang/project/xcited-state/qceims/MethylHistidine/GFnMethylHistidine.jdx'
p = [0.59, 0.41]
p = [i/100 for i in p]
y = readjdx(file)
dot0 = dot_sim(x,y)



#%%
state2 = path + 'result.jdx'
state1 = file
out = path + 'mix.jdx'
file = state2
y = readjdx(file)
dot1 = dot_sim(x,y)


#%%

y1 = readjdx(state1)
y2 = readjdx(state2)

y2[0:15]=0



ymix = p[0]*y1 + p[1]*y2
ymix = ymix/np.amax(ymix)*1000


dot_mix = dot_sim(x,ymix)

print('dot0, dot1, dot_mix')
print(dot0, dot1, dot_mix)
#%%
n = np.count_nonzero(ymix)
index = np.nonzero(ymix)[0]
intensity = ymix[index]

with open(out,'w') as f:
    f.writelines('##TITLE=Theoretical in-silico spectrum (MIX)'+'\n')
    f.writelines('##JCAMP-DX=4.24'+'\n')
    f.writelines('##DATA TYPE=MASS SPECTRUM')
    f.writelines('\n')
    f.writelines('##XUNITS=M/Z')
    f.writelines('\n')
    f.writelines('##YUNITS=RELATIVE INTENSITY')
    f.writelines('\n')
    f.writelines('##NPOINTS=' + str(n))
    f.writelines('\n')
    f.writelines('##PEAK TABLE=(XY..XY) 1')
    f.writelines('\n')
    for i in range(0,n):
        f.writelines(str(index[i]))
        f.writelines('    ' )
        f.writelines(str(round(intensity[i],2)))
        f.writelines('\n')
    f.writelines('##END=')
