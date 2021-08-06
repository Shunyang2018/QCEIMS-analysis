#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 15:20:12 2019

@author: shunyang
transform sdf, ie mainlib.sdf to csv file! Watch out they are very big!
"""
import pandas as pd
import numpy as np
import re

msp = '/Users/shunyang/project/TMS/TMS/TMS-NIST14/mainlib-tms.msp'
sdf = '/Users/shunyang/project/TMS/TMS/TMS-NIST14/mainlib-tms-filtered.sdf'
atoms = {'C', 'H', 'O', 'N', 'Si', 'S', 'P'}

molecule = 'C13H21NO2Si'
def mol_filter(molecule):

    
    
    All = set(re.findall(r'[A-Z][a-z]|[A-Z]',molecule))#here we need to get Si first, so we can't change the sequence of '|'
    return ('Si' in All) & (All.issubset(atoms))
 
mol_filter(molecule)


def read_blocks(file):
    block = []
    for line in file:
        if line.startswith('Name:') and len(block)>0:
            yield block
            block = []
            
        block.append(line)
    yield block
#i=1
#with open(msp) as f:
#    for block in read_blocks(f):
#        i += 1
#        if i == 4:
#            break


def mass2list(mass):
    x=np.zeros([700,])# mass from 0 to 600
    i = 0
    for line in mass:
        
        temp = line.strip().split(';')
        for m in temp[0:-1]:
            try:
                n = m.strip().split(' ')
#                print(n)
                i += 1
    #            print(int(temp[0]), float(temp[1]))
                x[int(n[0])] = float(n[1])
                
            except IndexError:
                break
            except ValueError:
                break
    return x,i
        
  

NAME = []
INCHIKEY = []
FORMULA = []
MASS = []
ExactMass = []
INCHI = []
timer = 1
#with open(sdf) as ff:
with open(msp) as f:
    for block in read_blocks(f):
#        print(timer,'compounds')
        timer += 1
        n = len(block)
        flag = False
        for i in range(0,n):
            if 'Name' in block[i]:      
                name = block[i].split(':')[1]
            if 'InChIKey' in block[i]:
                inchikey = block[i].split(':')[1]
                
            if 'Formula' in block[i]:
                print(block[i].split(':')[1])
                mol = block[i].split(':')[1]
                if mol_filter(mol):
                    flag = True
                    
            if 'Computed-InChI:' in block[i]:
                inchi = block[i].split(':')[1]
            if 'ExactMass' in block[i]:
                exact = block[i].split(':')[1]
            if ('Num Peaks' in block[i])&flag:
                Num = block[i].split(':')[1]
                mass = block[i+1:-1]
                test,N = mass2list(mass)
            if ('MW' in block[i]) & flag:
                mw = int(block[i].split(':')[1])
                if mw >= 700: # mass more than index
                    flag = False
        if flag == True:
            NAME.append(name)
            INCHIKEY.append(inchikey)
            ExactMass.append(exact)
            FORMULA.append(mol)
            if N != int(Num):
                
                print ('wrong!!!' + inchikey + str(N) + str(Num))
                break
            MASS.append(test)
            INCHI.append(inchi)
#                ff.write(block)
            
#INCHIKYEY = INCHIKEY.pop()           
print('NAME:',len(NAME),'INCHIKEY',len(INCHIKEY), 'FORMULA', len(FORMULA), 'MASS', len(MASS))               
#                
df = pd.DataFrame({'NAME':NAME, 'INCHIKEY':INCHIKEY, 'FORMULA':FORMULA,  'ExactMass':ExactMass, 'INCHI':INCHI   })
#df['SHORT'] = df['INCHIKEY'].str.split('-', n=1,expand=True)
df.to_csv('/Users/shunyang/project/TMS/TMS/TMS-NIST14/TMS_under700.csv')