#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 10:39:15 2021

@author: shunyang
"""

from decimal import Decimal
import os 


def D2float(string):
    tmp = string.split('D')
    ten = 10 ** int(tmp[1])
    return float(tmp[0])*ten

def float2D(digit):
    
    tmp = '%.14E' % Decimal(str(digit))
    tmp0 = tmp.split('E')[0]
    tmp1 = tmp.split('E')[-1]
    
    de = float(tmp0)/10
    ten = int(tmp1) + 1
    if digit > 0: 
        tmp = ' %f00000000D%+d'%(de,ten)
        
    else:
        tmp = '%f00000000D%+d'%(de,ten)
    return tmp.replace('D+','D+0').replace('D-','D-0')



def read_start(path):
    name = path + 'qceims.start'
    with open(name,'r') as f:
        firstline = f.readline()
        natom = int(firstline.split()[1])
        secondline = f.readline()
        tmp = secondline.split()
        iee = (D2float(tmp[0]))
        tadd = (D2float(tmp[1]))
        trunk = f.readlines()
        f.close()
        
    coord = []
    velo = []
    for i in range(len(trunk)):
        
          
        f=trunk[i].split()
        coord.append(list(map(D2float, f[0:3])))
        velo.append(list(map(D2float, f[3:6])))
    
    import numpy as np
    coord = np.array(coord)
    velo = np.array(velo)
    return coord, velo, [iee,tadd]
# def write_start(path):
def write_start(first_line, path, add, coord, velo):
    outfile = path + 'qceims.start'
    with open(outfile,'w') as f:
        f.write(first_line)
        secondline = (' '+float2D(add[0]) + '  '+ float2D(add[1])+ '\n')
        f.write(secondline)
        
        for n in range(len(coord)):
            string1 = ' '.join([float2D(i) for i in coord[n]])
            string2 = ' '.join([float2D(i) for i in velo[n]])
            f.write(' '+string1 + ' ' + string2 + '  0.10000000000000D+01\n')
            
#%% get atom list

tmol = '/Users/shunyang/project/xcited-state/wigner/cw_vn/coord'
name = []
with open(tmol, 'r') as r:
    tmp = r.readlines()
    for i in tmp[1:-1]:
        name.append(i.split()[-1].upper())               
#%% rmsd of two structures
from rmsd import kabsch_weighted_rmsd 

TMP = '347'
path = '/Users/shunyang/project/xcited-state/wigner/UA_wigner_py'
path1 = path + '/TMPQCEIMS/TMP.' +TMP+  '/'
path = '/Users/shunyang/project/xcited-state/wigner/non_wig'
path2 = path + '/TMPQCEIMS/TMP.' +TMP+ '/'


coord1, velo1,add1 = read_start(path1)
coord2, velo2,add2 = read_start(path2)
rmcoord = kabsch_weighted_rmsd(coord1,coord2)
rmvelo = kabsch_weighted_rmsd(velo1,velo2)
print(path1.split('.')[-1],
      'coord', 'velo')
print(' ',rmcoord, rmvelo)
print('path1',add1)
print('path2',add2)




first_line = ' '+TMP+'  16\n'
path = '/Users/shunyang/project/xcited-state/wigner/cw_vn/'
add = add1
coord = coord1.tolist()
velo = velo2.tolist()


write_start(first_line, path, add, coord, velo)


path = '/Users/shunyang/project/xcited-state/wigner/cn_vw/'
add = add1
coord = coord2.tolist()
velo = velo1.tolist()


write_start(first_line, path, add, coord, velo)    
    
#%% write xyz to validate the structures, especially about wigner distribution
def write_xyz(path, coord):
    coord = coord/1.8
    outfile = path + 'coord.xyz'
    with open(outfile,'w') as f:
        f.write('16 \n\n')
        
        for n in range(len(coord)):
            f.write('%-2s %22.15f %22.15f %22.15f\n' % ( name[n], coord[n][0], coord[n][1], coord[n][2]))
TMP = '184'
path = '/Users/shunyang/project/xcited-state/wigner/UA_wigner_py'
path1 = path + '/TMPQCEIMS/TMP.' +TMP+  '/'
coord1, velo1,add1 = read_start(path1)
write_xyz(path1, coord1)

