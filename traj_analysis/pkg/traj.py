#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 22:53:17 2019

@author: shunyang
"""

import numpy as np
#import tarfile
import re
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
rad = [       0.643,0.643,2.457,1.909,1.587,1.436,1.309,
       1.096,1.120,0.945,2.986,2.646,2.400,2.192,
       2.060,1.890,1.795,1.701,3.836,3.288,2.721,
       2.494,2.305,2.230,2.211,2.211,2.192,2.173,
       2.211,2.362,2.381,2.305,2.268,2.192,2.154,
       2.116,4.082,3.609,3.061,2.740,2.532,2.457,
       2.400,2.362,2.362,2.419,2.532,2.797,2.721,
       2.665,2.646,2.570,2.513,2.476,4.441,3.742,
       3.194,3.118,3.118,3.099,3.080,3.061,3.496,
       3.042,3.005,3.005,2.986,2.967,2.948,2.948,
       2.948,2.721,2.532,2.457,2.419,2.381,2.400,
       2.457,2.532,2.816,2.797,2.778,2.759,2.759,
       2.740,2.800,
       2.800,2.800,2.800,2.800,2.800,2.800,2.800,
       2.800,2.800,2.800,2.800,2.800,2.800,2.800]
rad = np.asarray(rad)
rcut = 3.0
at1 = 1
at2 = 0
elem = np.asarray(['h','he',
      'li','be','b','c','n','o','f','ne',
      'na','mg','al','si','p','s','cl','ar',
      'k','ca','sc','ti','v','cr','mn','fe','co','ni','cu',
      'zn','ga','ge','as','se','br','kr',
      'rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag',
      'cd','in','sn','sb','te','i','xe',
      'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
      'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
      'au','hg','tl','pb','bi','po','at','rn',
      'fr','ra','ac','th','pa','u','np','pu','am','cm','bk','cf','xx',
      'fm','md','cb','xx','xx','xx','xx','xx'])
np.where(elem == 'he')
#def getelem(name):
#    if name.lower() in elem:
#        return

def process(coord):

    cartesian = np.zeros((len(coord)-2,3))
    initial_line = float(re.split(r'\n',coord[1])[0]) # initial line contains several important info about the point
    atom_list = []

    for i in range(2,len(coord)):# extract out the coordinate of each atom and put it into a list

        line = re.split(r'[\s]',coord[i])

        tmp = re.findall(r'[A-Z][a-z]|[A-Z]',coord[i])[0].lower()

        atom_list.append(int(np.where(
                elem == tmp[0])[0])+1)

        line = [x for x in line if x!= '']
        cartesian[i-2,0] = float(line[1])
        cartesian[i-2,1] = float(line[2])
        cartesian[i-2,2] = float(line[3])

    return cartesian,initial_line,atom_list
def read_blocks(file):#iterate blocks with end sign"initial Cartesian coordinates:"
    block = []
    for line in file:

        if (not re.search(r'[*.*]', line)) and len(block)>0:
            yield block
            block = []
        block.append(line)
    yield block

def readxyz(file):
    tmptraj = []
    tmpenergy = []
    tmpatom = []
    with open(file) as f:
        for block in read_blocks(f):

            tmp = process(block)
            tmptraj.append(tmp[0])
            tmpenergy.append(tmp[1])
            tmpatom.append(tmp[2])

    f.close()
    return tmptraj,tmpenergy,tmpatom

def fragment(coord,iat):

    iat = np.array(iat)-1

    atom_num = coord.shape[0]

    connect = np.zeros((atom_num, atom_num))#0 is not connected
    r = np.zeros((atom_num, atom_num))
    frag= np.zeros(atom_num)
    for i in np.arange(0,atom_num-1,1):
        for j in np.arange(i+1,atom_num,1):
            r[i][j] = np.sqrt((coord[i][0]-coord[j][0])**2+
                        (coord[i][1]-coord[j][1])**2+
                        (coord[i][2]-coord[j][2])**2)
            r[j][i] = r[i][j]

            rcov = rcut*0.5*(rad[iat[i]] + rad[iat[j]])/0.52917
            if r[i][j] < rcov :
                connect[i][j] = 1
                connect[j][i] = 1

#    return connect, r

    if (at1 == 0)&(at2 == 0):
        for i in np.arange(0,atom_num,1):
            frag[i] = 1
        return frag
    else:
        for i in np.arange(0,atom_num,1):
            frag[i] = 0
        frag[at1-1] = 1
        attotal = 1
        if at2 != 0:
            connect[at1][at2] = 0
            connect[at2][at1] = 0
        flag = False
        currentfrag = 0



        while attotal != atom_num:

            currentfrag=currentfrag + 1

# cycle through atoms and find connected ones
            while not flag:
                flag = True
                for i in np.arange(0,atom_num-1,1):
                    if frag[i] == currentfrag:
                        for j in np.arange(i+1,atom_num,1):
                            if connect[i][j] == 1:
                                if frag[j] == 0:
                                    frag[j] = currentfrag
                                    attotal = attotal + 1

                                    flag = False
#                                elif frag[j] == currentfrag:
#                                    continue


            #find
            for i in np.arange(0,atom_num,1):
                if frag[i] == 0:
                    frag[i] = currentfrag + 1
                    attotal = attotal + 1
                    break
            flag = False

    return frag




class trajectory:
    def __init__(self,path):
        '''
        get trj.1.1 file name
        as a list
        '''
        self.path = path
        self.TMP_list = os.listdir(path)
        self.traj_namelist = []
        for i in self.TMP_list:
            i = path+i+'/'
            tmp = os.listdir(i)
            traj = []
            for j in tmp:
                if 'trj' in j:
                    traj.append(i+j)
            self.traj_namelist.append(traj)
    def out(self):
        return


#        tar_file = path
#        tar = tarfile.open(tar_file, "r:gz")
#        self.tar = tar.getmembers()
    def readtraj(self):
        '''
        save cartesian, energy and atom list
        '''
        self.coord = []
        self.energy = []
        self.iat = []
        for i in self.traj_namelist:
            coord1, energy1, iat1 = readxyz(i[0])#only look into the first traj, regardless further fragmentation

            self.coord.append( coord1)
            self.energy.append(energy1)
            self.iat.append(iat1)

        return
    def getfrag(self):
        self.frag = []
        self.nfrag = []
        for i in range(len(self.coord)):
            fragj = []
            nfragj = []
            for j in range(len(self.coord[i])):
                tmp = fragment(self.coord[i][j], self.iat[i][j])
                fragj.append(tmp)
                nfragj.append(np.max(tmp))
#                self.nfrag
            self.frag.append(fragj)
            self.nfrag.append(nfragj)
        return self.frag

    def plot(self):
        fig, ax = plt.subplots()
        for i in self.nfrag:
            ax.plot(2*np.arange(len(i)),i)
        ax.set_ylabel('Number of fragments')
        ax.set_xlabel('Time/fs')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # plt.show()
        plt.savefig('fragments.png')
        return
    def plotenergy(self):
        fig, ax = plt.subplots()
        for i in self.energy:
            ax.plot(2*np.arange(len(i)),i)
        ax.set_ylabel('Energy')
        ax.set_xlabel('Time/fs')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # plt.show()
        plt.savefig('energy.png')
        return
#%%
t = trajectory('/Users/shunyang/project/qceims_input/examples/newversion/184/TMPQCEIMS/')
t.readtraj()
test = t.coord[0]

test2 = t.getfrag()

test3 = t.nfrag
t.plot()
t.plotenergy()

test = t.energy
