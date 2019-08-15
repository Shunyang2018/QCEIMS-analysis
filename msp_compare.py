#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 15:31:58 2019

@author: shunyang
compare sigle msp file with our jdx results
"""


'''
run this section to load functions first
'''
#for NIST17
import pandas as pd
import numpy as np
import os


debug = True

  

def cos_sim(a, b):
    '''
    calculate cosine similarity of two vetcors
    Input
        a: np.array
        b: np.array
    Output
        res: float
                
    '''
    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    res = dot_product / (norm_a * norm_b)
    return res 


def dot_sim(x,y):
    '''
    calculate dot product of two mass spectrums
    Input 
        a: np.array
        b: np.array
    Output
        dot: float normalized in 1000
    '''
    m=0.6
    n=3
   
    #mast use different parameters!!!otherwise you will destory the initial array!!!!!!
    a = np.copy(x)
    b = np.copy(y)
    for i in range(0, len(a)-1):
        a[i] = (i**n)*(a[i]**m)
    for i in range(0, len(b)-1):
        b[i] = i**n*(b[i]**m)
    dot= np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))*1000
#    dot.append( np.power(np.dot(a,b)**2/(np.sum(a**2)*np.sum(b**2)), 0.5) )
    return dot  
    
    
    
    
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
    # another way to drop last column             
    #df = df.drop(df.columns[len(df.columns)-1],axis = 1)
    
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
    

def match(x, y):
    '''
    NIST match scores
    '''
    dot = dot_sim(x,y)
    N_y = np.count_nonzero(y)
    N = 0
    SUM = 0
    a = []
    b = []
    for i in range(1, len(y)):
        if (x[i] != 0) & (y[i] != 0):
            N += 1
            a.append(x[i])
            b.append(y[i])
    for j in range(1, N):
        temp = a[j]*b[j-1]/(a[j-1]*b[j])
        if temp > 1:
            temp = 1/temp
        SUM += temp
    print(SUM)
    return (N_y*dot + SUM*1000)/(N_y + N)


msp = '/Users/shunyang/Box/qceims_input/similarity/admantane695.MSP'#the experimental data
x = readmsp(msp)
path='/Users/shunyang/Box/qceims_input/conformer-2-nonene/adamantane_same' # .jdx file path 

#%%
'''
compare single jdx with msp
'''
file = '/Users/shunyang/Box/qceims_input/results/184.jdx'
y = readjdx(file)
y = y/np.amax(y)
dot = dot_sim(x,y)
cos = cos_sim(x,y)
N_y = np.count_nonzero(y)
N = 0
SUM = 0
a = []
b = []
m = 0
for i in range(1, len(y)):
    if (x[i] != 0) & (y[i] != 0):
        N += 1
        a.append(x[i])
        b.append(y[i])
for j in range(1, N):
    print(b[j-1])
    temp = a[j]*b[j-1]/(a[j-1]*b[j])
    if temp > 1:
        temp = 1/temp
    m += 1
    SUM += temp# pay attention to intent of for loop!!!

score = (N_y*dot + SUM*1000)/(N_y + N)

#%%
'''
jdx files in directory
output csv file of three kinds of similarity at the same path
'''
  
os.chdir(path)
files= os.listdir(path)
Dot = []
Cos = []
Match = []
for file in files:
    print(file)

    if os.path.splitext(file)[1] == ".jdx":
        y = readjdx(file)#read computational results
        y = y/np.amax(y)#normalization to 100
        Match.append(match(x,y))
#     here x is library and y is unknown
        cos = cos_sim(x,y)*1000
        dot = dot_sim(x,y)  
        Cos.append(cos)
        Dot.append(dot)    
    
dataframe = pd.DataFrame({'Index':files, 'Dot':Dot, 'Cos':Cos, 'Match':Match})
dataframe.to_csv('result.csv',index=False,sep=',')
#    
#%%
##### multiple MSP with same jdx
'''
compare one simulated result with multiple experimental results
'''
path = '/Users/shunyang/Box/qceims_input/similarity' 
jdx =  '/Users/shunyang/Box/qceims_input/results/184.jdx' 
x = readjdx(jdx) 
os.chdir(path)
files= os.listdir(path)
Dot = []
Cos = []
Match = []
for file in files:
    print(file)

    if os.path.splitext(file)[1] == ".MSP":
        y = readmsp(file)#read computational results
        y = y/np.amax(y)#normalization to 100
        Match.append(match(x,y))
#     here x is library and y is unknown
        cos = cos_sim(x,y)*1000
        dot = dot_sim(x,y)  
        Cos.append(cos)
        Dot.append(dot)       
    
    
    
    
    
    
