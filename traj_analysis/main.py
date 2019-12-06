#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 20:22:09 2019

@author: shunyang
"""
import pkg.traj as traj

t = traj.trajectory('/Users/shunyang/project/qceims_improve/traj_analysis/184/TMPQCEIMS/')
t.readtraj()
test = t.coord
test2 = t.getfrag()
