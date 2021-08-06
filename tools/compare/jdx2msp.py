#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 15:55:14 2021

@author: shunyang
"""

import os
import pandas as pd
from jdxtool import cos_sim, dot_sim, readjdx, peak, match,  trans
import numpy as np

jdxpath = '/Users/shunyang/project/coauthor/csv/msp/'
os.chdir(jdxpath)
files= os.listdir(jdxpath)
indexfile = '/Users/shunyang/project/TMS/TMS/TMS-NIST14/TMS_under700_classify_v5_correctmass.xlsx'
sheet = pd.read_excel(indexfile,sheet_name='single derivative reference', engine='openpyxl',index_col=None, usecols = "A:P",header=0)
NIST = '/Users/shunyang/project/TMS/TMS/TMS-NIST14/mainlib-tms.msp'



test = trans(readjdx('/Users/shunyang/project/TMS/TMS/spectra/GFn200.jdx'))
output = '/Users/shunyang/project/TMS/jdx2.msp'


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
    mass = df.values[:,0].astype(float).astype(np.int32)
    intensity = df.values[:,1].astype(np.int32)
    #generate computational spec
    x=np.zeros([699,])# mass from 0 to 699
    for i in range(0,len(mass)-1):
        x[mass[i]] = intensity[i]
    print('max',np.amax(x))
    x = x/np.amax(x)

    return x




accurate = False


with open(NIST,'r') as f:
    with open(output,'w') as w:
        n=0
        for file in files:
            # read jdx
            if 'jdx' in os.path.splitext(file)[1]:
                filename = file
                df = pd.read_csv(filename, header=6)
                different_list = df.iloc[:-1,0].str.split().str
                if not accurate:
                    df['Mass'],df['intensity']= df.iloc[:-1,0].str.split().str
                    df['formula'] = ' '
                    file = file.replace('result','')
                else:
                    file = file.replace('accurate','')
                    df['Mass'],df['intensity'],df['formula']= df.iloc[:-1,0].str.split().str
                df = df.drop(columns='##PEAK TABLE=(XY..XY) 1').drop(df.tail(1).index)
                file = file.replace('accurate','')
                # df = df.drop(['formula'],axis=1)
                Num = 'Num Peaks: ' + str(len(df))


                if os.path.splitext(file)[1] == '.jdx': # find jdx file
                    msindex = os.path.splitext(file)[0].replace('GFn','')
                    msindex = int(msindex)


                    inchikey = sheet['INCHIKEY'].loc[sheet['No.']==msindex].item() # read inchikey from excel file
                    print('working on...')
                    print(file,inchikey)



                    # read head information from NIST 17
                    data = f.read()
                    i = data.index(inchikey)
                    f.seek(i,0)
                    print('inchikey validate as:',f.readline())
                    i1=data.index('Num Peaks:',i,-1)
                    f.seek(i1,0)

                    i0 = data.rfind('Name:',0,i)
                    f.seek(i0,0)
                    data = data[i0:i1].splitlines()
                    # clean the head
                    for line in data:
                        if 'DB#' in line:
                            print(line)
                            data.remove(line)
                            break
                    DB = 'DB#: QCEIMS' + str(msindex)
                    data.insert(-1,DB)
                    other = 'Instrument_type: in-silico QTOF\nIon_mode: P\nSpectrum_type: MS1'
                    data.insert(-1,other)
                    print('writing file...')
                    for line in data:
                        if 'Comments:' not in line:
                            w.write(line+'\n')

                    w.write(Num+'\n')
                    w.write(df.to_string(header=False,index=False, index_names=False))
                    w.write('\n\n')
                    f.seek(0)
                    n += 1
                    print('n',n)