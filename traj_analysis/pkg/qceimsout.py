#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:55:50 2019

@author: shunyang
"""

import pandas as pd

'''
To use the code, change the path of qceims.out file, 
the program will generate a qceims.csv file at the same path
'''

path = '/Users/shunyang/project/qceims_input/184/qceims.out'

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
                self.time.append(block[i-2].split()[1])
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
                               'step   time [fs]    Epot       Ekin       Etot    error  #F   eTemp   frag. T': self.content}
        
        self.pd = pd.DataFrame(data)
        self.pd.to_csv(path.split('.')[0]+'.csv')
        return self.pd
                    
            


y = qceimsout(path)
print(len(y.trajectory1),len(y.trajectory2),len(y.exit),len(y.content))

y.excel()










