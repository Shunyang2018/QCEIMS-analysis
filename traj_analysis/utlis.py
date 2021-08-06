#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 11:09:58 2021

@author: shunyang
"""
import os 
import pandas as pd 
from collections import Counter
symbol = {'1':'H', '6':'C', '7':'N','8':'O','14':'Si'
          ,'9':'F', '16':'S','17':'Cl'}

res = 'qceims.res'

res_list = []
with open(res,'r') as f:
    res = f.readlines()
    for line in res:
        line = line.split()
        intensity = line[0]
        trj = line[1]
        frag0 = line[2]
        frag1 = line[3]
        number = int(line[4])
        mol = []
        for i in range(0,number):
            atom = [symbol.get(line[5+2*i]),line[6+2*i]]
           
            mol.append(atom)
            # break
        mol = sorted(mol, key=lambda x: (x[0]))
        formula = ''
        for pair in mol:
            formula += pair[0]+pair[1]
        res_list.append([intensity, trj, frag0, frag1, formula, mol])
        
    
    
    

# def parse_molecule(molecule):

#     array = [[]]
    
#     for token in map(str, molecule):
#         if token.isalpha() and token.istitle():
#             last = [token]
#             upper = token
#             array[-1].append(token)
#         elif token.isalpha():
#             last = upper + token
#             array[-1] = [last]
#         elif token.isdigit():
#             array[-1].extend(last*(int(token)-1))
#         elif token == '(' or token == '[':
#             array.append([])
#         elif token == ')' or token == ']':
#             last = array.pop()
#             array[-1].extend(last)

#     return dict(Counter(array[-1]))



import re
def parse_molecule(molecule):
    return dict(re.findall(r'([A-Z][a-z]?)(\d*)', molecule))



test = "Si40C2O"
out = parse_molecule(test)
out.get('O')

#%%

import openbabel


obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("xyz", "pdb")

trajectory = '/Volumes/File_Backup/backup/TMS/TMS-leucine/305/TMPQCEIMS/TMP.10/trj.10.2'
mol = openbabel.OBMol()
obConversion.ReadFile(mol, trajectory)
obConversion.WriteFile(mol, trajectory+'_pdb')
command = 'obabel -ixyz ' + trajectory + ' -opdb -O' + trajectory +'.pdb'
os.system(command)


#%%
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

import pytraj as pt
traj = pt.load('/Users/shunyang/project/TMS/Code/traj_analysis/trj.pdb')
data = pt.pca(traj, mask='!@H=', n_vecs=2)

print('projection values of each frame to first mode = {} \n'.format(data[0][0]))
print('projection values of each frame to second mode = {} \n'.format(data[0][1]))
print('eigvenvalues of first two modes', data[1][0])
print("")
print('eigvenvectors of first two modes: \n', data[1][1])

import matplotlib

projection_data = data[0]
from matplotlib import pyplot as plt

plt.scatter(projection_data[0], projection_data[1], marker='o', c=range(traj.n_frames), alpha=0.5)
plt.xlabel('PC1')
plt.ylabel('PC2')
cbar = plt.colorbar()
cbar.set_label('frame #')

df = pt.multidihedral(traj, dtype='dataframe')


#%%

import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import nglview as nv
import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')





u = mda.Universe('/Volumes/File_Backup/backup/TMS/TMS-leucine/305/TMPQCEIMS/TMP.10/trj.10.2.pdb')

aligner = align.AlignTraj(u, u,
                          in_memory=True).run()

pc = pca.PCA(u,
             align=False, mean=None,
             n_components=None).run()


backbone = u.select_atoms('backbone')
n_bb = len(backbone)
print('There are {} backbone atoms in the analysis'.format(n_bb))
print(pc.p_components.shape)

plt.plot(pc.cumulated_variance[:10])
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance');

transformed = pc.transform(u, n_components=5)
df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(5)])
df['Time (ps)'] = df.index * u.trajectory.dt

import seaborn as sns

g = sns.PairGrid(df, hue='Time (ps)',
                 palette=sns.color_palette('Oranges_d',
                                           n_colors=len(df)))
g.map(plt.scatter, marker='.')
