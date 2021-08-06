#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:10:10 2021

@author: shunyang
"""
# use pytraj for reading coordinates from rst7 file
import pytraj as pt

# use nglview for visualizing this notebook
import nglview as nv
from time import time, sleep

from IPython.display import display
import os
import openbabel
import re


class molecule():
    def __init__(self):
        pass
    def io(self, formula, path, mass):
        self.formula = self.parse_molecule(formula)
        self.path = path
        self.mass = int(mass)
        self.res = self.getres()

    def parse_molecule(self, formula):
        tmp = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        tmp = sorted(tmp,key=lambda x: (x[0]))
        tmp_clean = []
        for atom in tmp:
            if atom[1] == '':
                atom = (atom[0],'1')
            tmp_clean.append(atom)
        return tmp_clean
    def getres(self):
        symbol = {'1':'H', '6':'C', '7':'N','8':'O','14':'Si'
          ,'9':'F', '16':'S','17':'Cl'}
        atommass = {'H':1, 'C':12, 'N':14,'O':16,'Si':28
          ,'F':19, 'S':32,'Cl':35}
        res = self.path + '/qceims.res'

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
                    atom = (symbol.get(line[5+2*i]),line[6+2*i])

                    mol.append(atom)

                mol = sorted(mol, key=lambda x: (x[0]))
                formula = ''
                mass = 0
                for pair in mol:
                    formula += pair[0]+pair[1]
                    mass += atommass.get(pair[0]) * int(pair[1])
                 

                res_list.append([intensity, trj, frag0, frag1, formula, mol, mass])
        return res_list
    def search(self):
        found = []
        for entry in self.res:
            if self.formula in entry:
                found.append(entry[0:5])
        if len(found) != 0:
            print('%i trajecoties are found'%(len(found)))
            print('intensity  trj  frag0  frag1')
            self.found = found
            return found
        else:
            print(self.formula)
            print('formula not found in res file!')
    def search_mass(self):
        found = []
        for entry in self.res:
            if self.mass in entry:
                found.append(entry[0:5])
        if len(found) != 0:
            print('%i trajecoties are found'%(len(found)))
            print('intensity  trj  frag0  frag1')
            return found
        else:
            print(self.mass)
            print('mass not found in res file!')
    def choosetrj(self, trj):
        display('please choose trajectory you want to see')
        self.trj = trj



    def viewtrj(self):
        path = self.path+'/TMPQCEIMS/TMP.'+self.trj

        file_list = os.listdir(path)
        probes = []

        for file in file_list:
            if ('trj' in file)&('pdb' not in file):

                trajectory = path+'/'+file
                command = 'sed -i \'\' \'s/SI/Si/g\' ' + trajectory
                print('running comman:',command)
                os.system(command)
                command = 'obabel -ixyz ' + trajectory + ' -opdb -O' + trajectory +'.pdb'
                print('running comman:',command)
                os.system(command)
                #mol = openbabel.OBMol()
                #obConversion.ReadFile(mol, trajectory)
                #obConversion.WriteFile(mol, trajectory+'.pdb')

                probes.append(trajectory+'.pdb')
                #traj = pt.load(trajectory+'.pdb')
                #self.view.add_trajectory(traj)
        self.view = []
        for p in probes:
            traj = pt.load(p)
            self.view.append(nv.show_pytraj(traj))
            #self.view.add_trajectory(traj)

    def get_IEE(self):
        iee = []
        tadd = []
        for i in self.found:
            TMP  = self.path + '/' + 'TMPQCEIMS/TMP.'
            i = TMP + i[1] + '/qceims.start'
            with open(i,'r') as f:
                f.readline()
                tmp = f.readline().split()
                iee.append(tmp[0])
                tadd.append(tmp[1])
                f.close()
        return iee, tadd
    def get_component_index(self, component_name):
        """
        Given the name of the loaded file (or a partial string),
        returns the component index if loaded
        """
        all_components = self.view._ngl_component_names
        component_idx = all_components.index([i for i in all_components
                                              if component_name in i][0])
        return component_idx

    def threshold_isosurface(self, probe_name, threshold):
        """
        Controls the representation of the map for specific probe
        :param probe_name: string, one of 'donor', 'acceptor', 'apolar'
        :param threshold: threshold value for displaying the hotspot map
        """

        comp_idx = self.get_component_index(probe_name)
        repr_params = [{'type': 'surface',
                        'params': {'opacity': 0.4,
                                   'isolevelType': 'value',
                                   'isolevel': threshold,
                                   }}]
        self.view.set_representations(repr_params, component=comp_idx)