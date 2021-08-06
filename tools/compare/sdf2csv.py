#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:20:12 2019

@author: shunyang
transform sdf, ie mainlib.sdf to csv file! Watch out they are very big!
"""
import pandas as pd
import numpy as np


sdf = '/Users/shunyang/project/TMS/TMS/TMS-NIST14/mainlib-tmssdf.msp'


def read_blocks(file):
    block = []
    for line in file:
        if line.startswith('$$') and len(block)>0:
            yield block
            block = []
        block.append(line)
    yield block

def mass2list(mass):
    x=np.zeros([700,])# mass from 0 to 500

    for line in mass:
        line = line.strip('\n')
        line = line.split(';')
        for pair in line:
            pair = pair.split()
            if len(pair) == 2:
                if int(pair[0])>=700:
                    x = 'Out of range'
                    break
                else:
                    x[int(pair[0])] = int(pair[1])
    return x


def readsdf(sdf):
    NAME = []
    INCHIKEY = []
    FORMULA = []
    MASS = []
    timer = 1
    with open(sdf) as f:
        for block in read_blocks(f):
    #        print(timer,'compounds')
            timer += 1
            print('reading',timer)
            n = len(block)
            flag = False
            for i in range(0,n):
                if 'Name' in block[i]:
                    print(block[i])
                    NAME.append(block[i].split(':')[1])
                if 'InChIKey' in block[i]:
                    INCHIKEY.append(block[i].split(':')[1])
                    flag = True
                if 'Formula:' in block[i]:
                    FORMULA.append(block[i].split(':')[1])
                if 'Num Peaks:' in block[i]:
                    mass = block[i+1:]
                    test = mass2list(mass)
                    MASS.append(test)
            if flag == False:
                print(timer, 'No-Inchikey')
                INCHIKEY.append('No-Inchikey')
    # print()            
    INCHIKYEY = INCHIKEY.pop()
    # df = pd.DataFrame({'NAME':NAME, 'INCHIKEY':INCHIKEY, 'FORMULA':FORMULA, 'MASS':MASS})
    return {'NAME':NAME, 'INCHIKEY':INCHIKEY, 'FORMULA':FORMULA, 'MASS':MASS}

# test = readsdf(sdf)  
# NAME = []
# INCHIKEY = []
# FORMULA = []
# MASS = []
# timer = 1
# with open(sdf) as f:
#     for block in read_blocks(f):
# #        print(timer,'compounds')
#         timer += 1
#         print('reading',timer)
#         n = len(block)
#         flag = False
#         for i in range(0,n):
#             if 'Name' in block[i]:
#                 print(block[i])
#                 NAME.append(block[i].split(':')[1])
#             if 'InChIKey' in block[i]:
#                 INCHIKEY.append(block[i].split(':')[1])
#                 flag = True
#             if 'Formula:' in block[i]:
#                 FORMULA.append(block[i].split(':')[1])
#             if 'Num Peaks:' in block[i]:
#                 mass = block[i+1:]
#                 test = mass2list(mass)
#                 MASS.append(test)
#         if flag == False:
#             print(timer, 'No-Inchikey')
#             INCHIKEY.append('No-Inchikey')

# INCHIKYEY = INCHIKEY.pop()
# print('NAME:',len(NAME),'INCHIKEY',len(INCHIKEY), 'FORMULA', len(FORMULA), 'MASS', len(MASS))

# df = pd.DataFrame({'NAME':NAME, 'INCHIKEY':INCHIKEY, 'FORMULA':FORMULA, 'MASS':MASS  })
# #df['SHORT'] = df['INCHIKEY'].str.split('-', n=1,expand=True)
# df.to_csv('mainlibTMS.csv')