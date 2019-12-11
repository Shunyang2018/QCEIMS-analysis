#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:47:28 2019

@author: shunyang
"""

import argparse
import os


parser = argparse.ArgumentParser(prog='qceimsout',
                                 description='Parser of qceims.out')
parser.add_argument("-f",'--file',dest='file', action='store',
                    help='folder path of msp files, example: path = \'/tmp/msp\'')
parser.add_argument("-o",'--output',dest='output',action='store',
                    help='path of output file, default is the upper folder of \
                    input file, example: path = \'/tmp/merge.msp\'')
parser.add_argument('-e', '--eof',dest='end', action='store',default='\n',
                    help='the marker between each msp, the default setting is ENTER ')
args = parser.parse_args( )


#['-f','/Users/shunyang/project/qceims_input/msp' ]
#path = '/Users/shunyang/project/qceims_input/msp'
#output = '/Users/shunyang/project/qceims_input/merge.msp'
#end = '\n'
file = os.listdir(args.file)
if not args.output:
    args.output = args.file + 'merge.msp'
print('will read those files: \n', file)
os.chdir(args.file)
with open(args.output, 'a') as o:
    
    for msp in file:
        if '.MSP' in msp:
            print('reading %s'% msp)
            with open(msp) as f: 
#                block = 
                o.writelines(f.readlines())
                o.write(args.end)
                f.close()
                
    print('file recorded')        
    o.close()
                