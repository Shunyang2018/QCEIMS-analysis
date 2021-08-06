#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 13:36:23 2019
@author: shunyang
calculate rmsd between multiple files
To use, 'pip install rmsd' first.
For more information:
https://github.com/charnley/rmsd

"""



import os
import numpy as np
import pandas as pd
import argparse
#parser settings
debug=True
parser = argparse.ArgumentParser(prog='trajrmsd',
                                 description='calculate rmsd between multiple files')
parser.add_argument("-t",'--traj',dest='traj', action='store',
                    help='trajectory file of xyz coordinates collection , example: path = \'/tmp/184traj.xyz\'')
parser.add_argument("-b",'--batch',dest='path',action='store',
                    help='path of folder which contains all structure files, example: \
                    /tmp/conformer-2-nonene/GMMX/1.xyz, than path = \'/tmp/conformer-2-nonene/GMMX/\', careful about the structure!')
parser.add_argument("-r", '--reference',dest='ref',action='store',
                    help='full path of the reference structure xyz files')
parser.add_argument("-o", '--output',dest='out',action='store',
                    help='output file name')

if debug:
    args = parser.parse_args(['-b','/Users/shunyang/project/TMS/TMS/TMS-traj/4-chlorobutan-1-ol/startconformer',
                            '-r','/Users/shunyang/project/TMS/TMS/TMS-traj/4-chlorobutan-1-ol/conformer_search'])
else:
    args = parser.parse_args()
#use this to test: ['-f','/Users/shunyang/project/qceims_input/184/qceims.out']

def rmsd(xyz,b):

    list_rmsd=[]
    list_rmsd_re=[]
    list_rmsd_h=[]
    for n in np.arange(0, len(xyz)):
        #print(xyz[n])
        a=xyz[n]
        print(a)
        #rmsd=float(os.popen('calculate_rmsd %s %s'%(a,b)).read().split('\\')[0])
        rmsd_re=float(os.popen('calculate_rmsd --reorder %s %s'%(a,b)).read().split('\\')[0])
        if debug:

            print(rmsd_re)
        #rmsd_h=float(os.popen('calculate_rmsd --reorder --no-hydrogen %s %s'%(a,b)).read().split('\\')[0])
        #list_rmsd.append(rmsd)
        list_rmsd_re.append(rmsd_re)
        #list_rmsd_h.append(rmsd_h)
    dataframe = pd.DataFrame({'rmsd_re':list_rmsd_re})#, 'rmsd_re':list_rmsd_re, 'rmsd_h':list_rmsd_h
    return dataframe


#%% main program
if (not args.path) & (not args.traj):# if input option wrong
    parser.error('Where is coordinate file!')
xyz=[]
if not args.traj:
    # compare xyz files in the path
    os.chdir(args.path)
    files= os.listdir(args.path)
    for file in files:
        if os.path.splitext(file)[1] == ".xyz":
            file = args.path+'/'+file
            xyz.append(file)
    xyz.sort()

elif not args.path:
    print('next time...I\'m lazy now. You can use openbable to convert trajectory file \
            into multiple xyz files and than use -b mode')

if not args.ref:
    print('using the first structure:%s  \n as a reference to calculate rmsd'%(xyz[0]))
    ref=xyz[0]
    dataframe = rmsd(xyz,ref)
else:
    xyzr = []
    os.chdir(args.ref)
    files= os.listdir(args.ref)
    for file in files:
        if os.path.splitext(file)[1] == ".xyz":
            xyzr.append(file)
    xyzr.sort()
    dataframe=pd.DataFrame()
    dfs = []
    for r in xyzr:#calculate rmsd with each reference
        print('using the reference structure: %s \n\t as a reference'%(r))
        df1 = rmsd(xyz,r)
        dfs.append(df1)
    dataframe = pd.concat(dfs,axis=1,keys=xyzr)
#%% save output
if not args.out:
    out='rmsd_traj.csv'
else:
    out = args.out
print('save output file as',out)
dataframe.to_csv(out,index=False,sep=',')




