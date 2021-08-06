#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 00:18:19 2021

@author: shunyang
"""
# use pytraj for reading coordinates from rst7 file
import pytraj as pt

# use nglview for visualizing this notebook
import nglview as nv
from time import time, sleep
from ipywidgets import widgets, interactive, VBox


traj = pt.load('trj.pdb')
view = nv.show_pytraj(traj)
VBox([view])

