#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:55:50 2019

@author: shunyang
"""

import pandas as pd
import argparse

'''
To use the code, change the path of qceims.out file, 
the program will generate a qceims.csv file at the same path
'''

file = '/Users/shunyang/project/qceims_input/184/qceims.out'

parser = argparse.ArgumentParser(prog='shunyang',
                                 description='Progdyn main program.')
parser.add_argument("-f",'--file',dest='file', action='store',
                    help='full path of the qceims.out file')
parser.add_argument("-b",'--batch',dest='path',action='store',
                    help='path of folder which contains all result file, example: /tmp/184/qceims.out, than path = /tmp/')
args = parser.parse_args(['-f',file])

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
#        self.calls = []
        with open(path) as f:
            for block in self.readblocks(f):
                self.block = block
                self.process(block)
                
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
#        flag= False
        for i in range(len(block)):
            if 'trajectory' in block[i]:
                tmp = block[i].split()
                self.trajectory1.append(tmp[1])
                self.trajectory2.append(tmp[2])
#                flag = not flag
            if ('E X I T' in block[i])|('EXIT' in block[i]):
                self.exit.append(block[i])
                self.content.append(block[i-2])
                tmp2 = block[i-2].split()
                self.time.append(tmp2[1])
                self.step.append(tmp2[0])
                self.Epot.append(tmp2[2])
                self.Ekin.append(tmp2[3])
                self.Etot.append(tmp2[4])
                self.error.append(tmp2[5])
                self.numfrag.append(tmp2[6])
                self.eTemp.append(tmp2[7])
                self.fragT.append(tmp2[8:])
            if 'M=' in block[i]:
                
                frag.append(block[i])
#                flag = not flag
#        if not flag:
#            self.re = block
#        if len(self.trajectory1) != len(self.exit):
#            print(block)
#            raise SystemExit('wrong!!!')
#            if 'QC calls' in block[i]:
#                self.calls.append(re.findall(r'*[0-9]',block[i]))
               
                
    def excel(self):
        data = {'index1': self.trajectory1, 'index2': self.trajectory2, 
                               'reason': self.exit, 'time': self.time,
                               'step':self.step, 'Epot':self.Epot, 'Ekin':self.Ekin,
                               'error':self.error, 'numfrag':self.numfrag, 'eTemp':self.eTemp, 'fragT':self.fragT}
        
        self.pd = pd.DataFrame(data)
        self.pd.to_csv(self.path.split('.')[0]+'.csv')
        return self.pd
                    
            


y = qceimsout(args.file)
print(len(y.trajectory1),len(y.trajectory2),len(y.exit),len(y.content))

y.excel()










