#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 14:46:43 2019

@author: shunyang
merge multiple mol file into one sdf file by the '$$$$' sign 
"""

import os
import sys



scriptPath = os.getcwd() + '/'
os.sys.path.append(scriptPath)

# open task file
task_filename = ('task-list.txt')# you need to prepare a list containing all the name of compounds
task_filename_handle = scriptPath + task_filename

# read all tasks from the file task-list.txt
my_tasks = ''

try:
    with open(task_filename_handle, "r") as myfile:
        my_tasks = myfile.readlines()
        myfile.close()
except IOError as e:
    print ('IO Error opening task-list.txt:  ', e)
    sys.exit(1)

# remove newline from list requires python 3
my_tasks = list(map(str.strip, my_tasks))
# printing the tasks
for x in range(len(my_tasks)):
    print (my_tasks[x])

# check if files can be written otherwise quit
try:
    with open(scriptPath + 'tmpfile.txt', "w") as myfile:
        myfile.write('')
        myfile.close()
except IOError as e:
    print ('IO Error writing file in directory:  ', e)
    sys.exit(1)


path = '/Users/shunyang/Box/qceims_input/database/'


with open('similarity_structure.sdf', 'w') as sdf:
    for name in my_tasks:
        filename = path + name + '.mol'
        with open(filename, 'r') as f:
            block = f.read()
            sdf.write(block)
            sdf.write('$$$$\n')
            f.close()            
    sdf.close()
    


