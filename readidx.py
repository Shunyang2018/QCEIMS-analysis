#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 15:31:58 2019

@author: shunyang
"""

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


def dot_sim(a,b):
    '''
    calculate dot product of two mass spectrums
    Input 
        a: np.array
        b: np.array
    Output
        dot: float normalized in 1000
    '''
    m=1.2
    n=0.9
    for i in range(0, len(a)-1):
        a[i] = (i**n)*(a[i]**m)
    for i in range(0, len(b)-1):
        b[i] = i**n*(b[i]**m)
    dot = 1000*(np.dot(a,b)**2)/(np.sum(a**2)*np.sum(b**2))
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
    x=np.zeros([501,])# mass from 0 to 500
    for i in range(0,len(mass)-1):
        x[mass[i]] = intensity[i]
    x = x/np.amax(x)            
    return x
    # another way to drop last column             
    #df = df.drop(df.columns[len(df.columns)-1],axis = 1)
def readdatabase(database):
    '''
    read database of experimental mass spectrum
    Input
        filename: str, name and path of the jdx file
    Output
        dataframe
    '''

    df = pd.read_excel(database,sheet_name=[0],header=0,index_col=0)[0]#mast have[0],to remove OrderedDict in header
    return df
    
    


def compare(filename,df):
    '''
    compare computational and experimental mass spectrums
    Input 
        filemane: str
            path and name of computational ms
        df: pandas.dataframe
            databases of experimental ms
    Output
        cos:float
            cos similarity
        dot:float
            dot product similarity like NIST search
    '''
    #get mass index in data base
    massindex = filename.split('/')[-1].split('.')[0]
    x = readjdx(filename)#read computational results

    if massindex not in df['Short'].tolist():
        print(massindex+'is not in database')
        cos='NA'
        dot='NA'
    else:
        spec=df[df['Short'] == massindex].index.tolist()[0] #get mass spectrum from massindex
#        index = df['mz 1'].index
#        print(index)
        y = np.array(df.loc[spec]['mz 1':].tolist())#from series to np.array, 19 is the first column of mass
        y = np.insert(y,0,0)#add 0 while mass=0, to get consistent with computational
        y = y/np.amax(y)#normalization to 100
        cos = cos_sim(x,y)
        dot = dot_sim(x,y)
    return cos, dot




database = '/Users/shunyang/Box/qceims_input/code/nist17.csv'    
filename = '/Users/shunyang/Box/analysising/MONA-EI-11240.jdx' 
output = '/Users/shunyang/Box/compared.csv'



    
path = "/Users/shunyang/Box/database" #file path
print("reading databases")
df = readdatabase(database)
print(df)
print("starting analysing")

#massindex = 'MONA-EI-11240'
#if massindex not in df['Short']:
#    print(massindex+'is not in database')
os.chdir(path)
files= os.listdir(path)
Dot = []
Cos = []
for file in files:
    print(file)
    cos, dot = compare(file,df)  
    Dot.append(dot)
    Cos.append(cos)
  
dataframe = pd.DataFrame({'Index':files, 'Dot':Dot, 'Cos':Cos})

if debug== True:

    print(dataframe)
dataframe.to_csv(output,index=False,sep=',')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    