#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:55:50 2019

@author: shunyang
anaylse qceims.out, give out running time, fragments, energy and some other informations
two modes: -f , for single file and -b for all out file in the secondary folder
"""

import pandas as pd
import argparse
import os 
import re

#file = '/Users/shunyang/project/qceims_input/184/qceims.out'

parser = argparse.ArgumentParser(prog='qceimsout',
                                 description='Parser of qceims.out')
parser.add_argument("-f",'--file',dest='file', action='store',
                    help='full path of the qceims.out file, example: path = \'/tmp/184/qceims.out\'')
parser.add_argument("-b",'--batch',dest='path',action='store',
                    help='path of folder which contains all result file, example: \
                    /tmp/184/qceims.out, than path = \'/tmp/\', careful about the structure!')
args = parser.parse_args( )
#['-f','/Users/shunyang/project/qceims_input/184/qceims.out']
    
#%%

def batchmode(path,batch):
    tmp = os.listdir(path)
    print('input list',tmp)
    for i in tmp:
        print('start to analyze: \n', i)
        tmppath=path + '/'+ i+'/qceims.out'
        print(tmppath)
        y = qceimsout(tmppath)

        ex = y.excel(batch)
        ex.to_csv(path + i + '.csv')
        print('output file at: \n', tmppath)

def count(coord):
    '''
    get the formula, could be wrong with Si!!!
    '''
    mol = ''
    for i in set(re.findall('[A-Z]', coord)):  
        
        mol += (i+str(coord.count(i)))
    return mol

def frag2list(x):
    '''
    clean the fragments column
    '''
    if type(x) != str:
        return ['','','','','','','','']
    else:
        x = x.replace('~','')
        x = x.replace('M=','')
        return x.split()         
#%%
class qceimsout():
    def __init__(self,path):
        '''
        get qceims.out file
        as a list
        '''
        self.path = path
        self.trajectory1 =[]
        self.trajectory2 = []
        self.exit = []
        self.content = []
        self.time = []
        self.step = []
        self.Epot = []
        self.Ekin = []
        self.Etot = []
        self.error = []
        self.numfrag = []
        self.eTemp = []
        self.fragT = []
        self.fraginfo = [] #only take the largest charge fragment
        self.frag = []
        self.start = []
        self.calls = []
        self.file(path)
        
#        else:
#            self.batch(path)
#        self.calls = []
    def file(self,path):
        '''
        iterator to generate a block including one TMP folder's qceims.out file
        '''
        print('reading files ...')
        with open(path) as f:
            for block in self.readblocks(f):
                self.block = block
                self.process(block)#read and record data
                if len(self.calls) != len(self.exit):
                    print('something wrong with data namber, \n has %i qcCalls'%len(self.calls))
#                    print(block)
                    break
#                break
                
        f.close()

            
    def readblocks(self, file):
        block = []
        flag = False
        for line in file:
            if 'ntraj' in line:
                flag = True
            if flag :
                if ('normal termination of QCEIMS' in line):
                    yield block
                    flag = False
                    block = []
                block.append(line)


    def process(self, block):
        frag = []
        flagcoord = False
        coord = ''
        timer = 0
#        flag= False
        for i in range(len(block)):
            if (flagcoord == True) & (len(block[i].split()) != 0):
#                print(block[i].split())
                coord += block[i].split()[4]
                
            elif (len(block[i].split())== 0) & (len(coord) != 0) :
                self.start.append(count(coord))
                coord = ''
                flagcoord = False
            if 'initial Cartesian coordinates:' in block[i]:
                flagcoord = True
            
            
                
             
            if 'trajectory' in block[i]:
                tmp = block[i].split()
                self.trajectory1.append(tmp[1])
                self.trajectory2.append(tmp[2])
                timer += 1
                
            if 'M=' in block[i]:
                
                frag.append(block[i].strip('\n'))
            if ('trajectory' in block[i]) & (len(frag) != 0):
#                print('frag',frag)
                self.frag.append(frag)
                frag = []
            if 'of QC calls' in block[i]:
                calls = block[i].split()[-1]
                self.frag.append(frag)
              
                
                for j in range(timer-1):
                    
                    self.calls.append(None)
                self.calls.append(calls)
                frag = []
                timer = 0 #don't know why we need!!! but with out this, the timer will be cumulated
#                flag = not flag
            if ('E X I T' in block[i])|('EXIT' in block[i]):
                self.exit.append(block[i].replace('E X I T   M D  because of ',''))
                self.content.append(block[i-2])
                if 'EXIT' in block[i]:
                    tmp2 = block[i-1].split()
                else:
                    tmp2 = block[i-2].split()
#                print(tmp2)
                self.time.append(tmp2[1])
                self.step.append(tmp2[0])
                self.Epot.append(tmp2[2])
                self.Ekin.append(tmp2[3])
                self.Etot.append(tmp2[4])
                self.error.append(tmp2[5])
                self.numfrag.append(tmp2[6])
                self.eTemp.append(tmp2[7])
                self.fragT.append(tmp2[8:])
            


    def list2column(self, oldname):
        '''
        seperate list in the column and rename them
        '''
        x = self.pd[oldname].apply(pd.Series)
        tmp = len(x.columns[:])
        name = []
        if len(oldname) > 9:
            x.rename(columns=dict(zip(x.columns[:],['mass', 'formula','index1','index2',
                                      'q pop','spin', 'q IDB', 'diss time'])),inplace=True)
            x = x.drop('index1',1)
            x = x.drop('index2',1)
        else:
                
            for i in range(tmp):
                name.append(oldname+str(i+1))
            x.rename(columns=dict(zip(x.columns[:],name)),inplace=True)
        
        self.pd = self.pd.drop(oldname , 1)
        self.pd = pd.concat([self.pd,x],axis=1)
        return self.pd              
              
    def excel(self,batch):
        data = {'index1': self.trajectory1, 'index2': self.trajectory2, 
                               'reason': self.exit, 'time': self.time,
                               'step':self.step, 'Epot':self.Epot, 'Ekin':self.Ekin,
                               'error':self.error, 'numfrag':self.numfrag, 
                               'eTemp':self.eTemp, 'fragT':self.fragT,'fragments':self.frag,'input':self.start,
                               'qcCalls':self.calls}
        
        self.pd = pd.DataFrame(data)
        self.list2column('fragT')
        self.list2column('fragments')
                
        for col in self.pd.columns:
            if 'fragments' in col:
             
                xx = self.pd[col].map(frag2list)
                self.pd = self.pd.drop(col , 1)
                self.pd = pd.concat([self.pd, xx],axis=1)
                self.list2column(col)
#        print(self.pd)
        if not batch:
            print('%i trajectories have been analyzed'%len(self.trajectory1))
            self.pd.to_csv(self.path.split('.')[0]+'.csv')
            print('output file at: \n', self.path.split('.')[0]+'.csv')
        return self.pd
                   
#%%           
if (not args.path) & (not args.file):
    parser.error('Where is qceims.out!')
if not args.path:
    batch = False
    print('start to analyze: \n', args.file)
    y = qceimsout(args.file)
    print(len(y.trajectory1),len(y.numfrag),len(y.exit),len(y.content),len(y.calls),
    len(y.time),len(y.step),len(y.Ekin),len(y.start),len(y.fragT))
    test = y.calls
#    print(y.frag)
    y.excel(batch)
    
elif not args.file:
    batch = True
    print('start batch mode!')
    batchmode(args.path,batch)


print('normal termination!')










