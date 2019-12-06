#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 20:24:15 2019

@author: shunyang
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


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
            x=np.zeros([499,])
            for line in block:
                line = line.strip('\n')
                line = line.split(';')
                for pair in line:
                    pair = pair.split()
                    if len(pair) == 2:
                        x[int(pair[0])] = int(pair[1])
    return x
def readjdx(filename):
    '''
    read and clean computational mass spectrum result jdx file
    Input
        filename: str, name and path of the jdx file
    Output
        x: normalized mass spectrum
    '''
    df = pd.read_csv(filename, header=6)
    
    df['Mass'],df['intensity']= df.iloc[:-1,0].str.split().str
    df = df.drop(columns='##PEAK TABLE=(XY..XY) 1').drop(df.tail(1).index)
    #get mass and intensity 
    mass = df.as_matrix()[:,0].astype(float).astype(np.int32)
    intensity = df.as_matrix()[:,1].astype(np.int32)
    #generate computational spec
    x=np.zeros([499,])# mass from 0 to 499
    for i in range(0,len(mass)-1):
        x[mass[i]] = intensity[i]

    x = x/np.amax(x) 
         
    return x

#Int = readmsp('/Users/shunyang/Box/qceims_input/similarity/admantane695.MSP')

def msplot(jdx):
    fig = plt.figure()
    Int = readjdx(jdx)
    name = jdx.split('/')[-1].replace('.jdx','')
    ax = fig.add_subplot(111)
    
    
    
    for i in range(len(Int)):
        if Int[i] != 0:
            h = int(Int[i]*1000)
    
            y=np.arange(0, h)
            x = np.ones(h)*i
            ax.plot(x, y, c ='b', label = name)
            ax.annotate(i,(i,h),xytext=(i, h+100 ),rotation=60,size=8,
#                        xycoords='axes points', rotation=90,
                        
#                        textcoords='axes points', size=8,
                        arrowprops=dict(arrowstyle="-",
                                        connectionstyle="Arc", linewidth=0.5, color='#808080'),
                        horizontalalignment='center', verticalalignment='bottom')
    #Get artists and labels for legend and chose which ones to display
    display = (0,499)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([handle for i,handle in enumerate(handles) if i in display],
               [label for i,label in enumerate(labels) if i in display])
    #ax.legend(name)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
#    ax.set_xlim(20,90)
    ax.set_ylim(0,1100)
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    plt.show()
    plt.savefig(name+'.png')
    
jdx = "/Users/shunyang/Box/qceims_input/parameters/2,4-di/143_tinit_300.jdx"
path = '//Users/shunyang/Box/qceims_input/conformer-2-nonene/adamantane/184conformer'
os.chdir(path)
files= os.listdir(path)

for file in files:
    if ('184' in file)&('jdx' in file):
            print (file)
            msplot(file)
            




#with open('admantane.xy','w') as f:
#    for i in range(len(Int)):
#        f.write(str(i)+' '+str(Int[i])+'\n')
#    f.close()
#        