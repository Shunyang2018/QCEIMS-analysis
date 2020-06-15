#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 08:42:10 2020

@author: shunyang

classify 1-TMS compounds by smiles, can be used through command line.


"""


from rdkit import Chem
from rdkit.Chem import AllChem, rdchem
from rdkit.Chem.rdchem import Atom 
import pandas as pd

def getAtom(m,element):
    '''
    get atom connected to the specific atom, which is 
    decided by GetIdx() or 'Si' elements.

    Parameters
    ----------
    m : rdkit.Chem.rdchem.Mol 
        the mol class containing the structure information
    element : str or int
        str only supports 'Si' now, and int is the atom index in m

    Returns
    -------
    TYPE list
    list of atom type and index connected to the specific atom,
    for example, (14, 'Si', [(13, 'O'), (15, 'C'), (16, 'C'), (17, 'C')])
    '''
    
    if element == "Si":
        # locate the atoms connected to Si
        for atom in m.GetAtoms():
        
            if atom.GetSymbol() == element:
                return atom.GetIdx(), atom.GetSymbol(), [(nbr.GetIdx(), nbr.GetSymbol()) 
                                for nbr in atom.GetNeighbors()]
    elif (type(element) == int):
        atom = m.GetAtomWithIdx(element)
#        print(atom.GetSymbol(),type(atom.GetIsAromatic()))
        return atom.GetIdx(), atom.GetSymbol(), [(nbr.GetIdx(), nbr.GetSymbol()) 
                            for nbr in atom.GetNeighbors()]

def Carbonyl(m,element):
    '''
    whether it is a C=O bond.

    Parameters
    ----------
    m : rdkit.Chem.rdchem.Mol 
        the mol class containing the structure information
    element : str or int
        str only supports 'Si' now, and int is the atom index in m

    Returns
    -------
    bool type

    ''' 
    Hybrid = mol.GetAtomWithIdx(element).GetHybridization()
    if Hybrid != rdchem.HybridizationType.SP2:
        
        return False
    else:
        tmp = getAtom(m,element)[2]
        # if len(tmp) > 2:
        #     print(tmp)
        for pair in tmp: #connected to alpha position atom 
            
            if pair[1] == "O":
                
                index2 = pair[0]
                
                bondtype = m.GetBondBetweenAtoms(element,index2).GetBondType()
                if bondtype == rdchem.BondType.DOUBLE:
                        return True #once find C=O return
        return False #no O atom at all

#test
# smi="CC(C)CC(C(=O)O[Si](C)(C)C)N"         
# mol = Chem.MolFromSmiles(smi)
# tmp = getAtom(mol,7)[2] 
# test = Carbonyl(mol,5)      

def betacarbons(mol, betalist):
    '''
    given molecule and beta atom list, judge whether containing C=O and in the aromatic ring

    Parameters
    ----------
    mol : TYPE
        DESCRIPTION.
    betalist : TYPE
        DESCRIPTION.

    Returns
    -------
    iscarbonyl : TYPE
        DESCRIPTION.
    isaromatic : TYPE
        DESCRIPTION.

    '''
    iscarbonyl = False # carbonyl flag
    isaromatic = False
    for pair in betalist: # go through all carbons, regardless carbon classification
        
        if pair[1] != "Si":
            if debug:
                print(pair)
            index2 = pair[0]
            atom2 = mol.GetAtomWithIdx(index2) #beta atom
            
            if atom2.GetIsAromatic(): 
                if debug:
                    print('aromatic',pair)
                isaromatic = True
            if Carbonyl(mol, index2):
                iscarbonyl = True
    return iscarbonyl, isaromatic
        

    
 #%%  parser part                
import argparse
debug=False
parser = argparse.ArgumentParser(prog='classification',
                                 description='classify 1-TMS compounds (must contains 1 TMS group), supporting alcohol, carboxylic acids, amines, amides, thiols now')
parser.add_argument('-b',dest='batch', action='store',help='path of smiles code conlumn, for example structure.smi')
parser.add_argument('-s',dest='smi', action='store',help='smiles string of given compounds, and will print classifications right way')
parser.add_argument('-o',dest='output',action='store',help='only activated for batch mode, save a csv file of classifications')
if debug:
    args = parser.parse_args( ['-b','/Users/shunyang/project/TMS/Functional_group/structure.smi'
                               #'-s','COc1ccccc1CCC(=O)O[Si](C)(C)C	' # primary acid
                               #'-s', 'Cn1c(cc2ccccc12)C(=O)O[Si](C)(C)C' #aromitic acid
                               #'-s', 'C/C(=N/OC)/[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4C[C@@H](CC[C@]4(C)[C@H]3CC[C@]12C)O[Si](C)(C)C' #secondary alcohol
                               #'-s', 'CC1(CCCCC1)O[Si](C)(C)C' #ter alcohol
                               #'-s','C/C(=N\O[Si](C)(C)C)/c1ccc(cc1)OC	'
                               ,'-o','/Users/shunyang/project/TMS/Functional_group/classifications.csv'])
else:
    args = parser.parse_args( )
print('program starting')


#%%
if args.batch:
    
    smifile= args.batch
    f = open(smifile, 'r')
    suppl = f.readlines()
elif args.smi:
    suppl = [args.smi]
if (not args.batch) & (not args.smi):# if input option wrong
    parser.error('missing parameters, use -h for help')
Si_list = []
type_list = []
subtype_list = []
test = []
for smi in suppl:
    mol = Chem.MolFromSmiles(smi)
    tmp = getAtom(mol,"Si")[2] # get the atoms connected to Si
    for xy in tmp:
        if xy[1] != "C": # get the only non-carbon atom on TMS-Si 
            atom = xy[1]
            index = int(xy[0])
    tmp2 = getAtom(mol,index) # beta position
    
    Si_list.append(tmp)
    # type_list.append(atom)
    atom = mol.GetAtomWithIdx(index)
    element = atom.GetSymbol()
    if element == "N":
        
        iscarbonyl, isaromatic = betacarbons(mol, tmp2[2])
        if iscarbonyl:
            type_list.append('amides')
        else:
            type_list.append('amines')
        if isaromatic:
            subtype_list.append("aromatic")
        else:
            if atom.GetIsAromatic():
                subtype_list.append("in aromatic ring")
            elif len(tmp2[2]) == 2:
                subtype_list.append("primary")
            elif len(tmp2[2]) == 3:
                subtype_list.append("secondary")
            else:
                subtype_list.append("unknown")
                print('unknown type',tmp2)

        
    elif element == "S":
        type_list.append('thiols')
        iscarbonyl, isaromatic = betacarbons(mol, tmp2[2])
        if isaromatic:
            subtype_list.append("aromatic") 
        else:
            if atom.GetIsAromatic():
                subtype_list.append("in aromatic ring")
            elif len(tmp2[2]) == 2:
                subtype_list.append("primary")
            elif len(tmp2[2]) == 3:
                subtype_list.append("secondary")
            elif len(tmp2[2]) == 4:
                subtype_list.append("tertiary")
            
    elif element == "O":
        iscarbonyl, isaromatic = betacarbons(mol, tmp2[2])
        tmp3=False #flag be flase if there is no Carbon on beta position, meaning O=N bond or someother kind.
        for xy in tmp2[2]: # all atom in beta position
            if xy[1] == "C": #to get  gamma carbon
                index2 = int(xy[0])
                tmp3 = getAtom(mol,index2)
        if tmp3:
            if iscarbonyl:
                type_list.append('acid')
                
                isaromatic2 = betacarbons(mol, tmp3[2])[1]
                if isaromatic2:
                    subtype_list.append("aromatic")
                else:
                    subtype_list.append("primary")
            else:
                type_list.append('alcohol')
                if isaromatic:
                    subtype_list.append("aromatic")
                else:
                    if atom.GetIsAromatic():
                        subtype_list.append("in aromatic ring")
                    elif len(tmp3[2]) == 2:
                        subtype_list.append("primary")
                    elif len(tmp3[2]) == 3:
                        subtype_list.append("secondary")
                    elif len(tmp3[2]) == 4:
                        subtype_list.append("tertiary")
            
        else:
            type_list.append('alcohol')
            subtype_list.append(tmp2[2])
    if debug:
        df = pd.DataFrame({ 'superclass': type_list, 'subclass' : subtype_list}) 
                
                # for pair in tmp2[2]: #connected to alpha position atom; beta position
                #     if pair[1] != "Si":
                #         index2 = pair[0]
                #         atom2 = pair[1]
                #         tmp3 = getAtom(mol,index2)[2] #beta position 
                        
                #         Hybrid = mol.GetAtomWithIdx(index2).GetHybridization()
                #         if Hybrid == rdchem.HybridizationType.SP3:
                #             subtype_list.append("alcohol")
                #         elif Hybrid == rdchem.HybridizationType.SP2:
                #             subtype_list.append("acid")
                #         else:
                #             subtype_list.append('unknown O')
                        
                        

        

        
        
#%%output part
        
        


df = pd.DataFrame({'smi': suppl, 'superclass': type_list, 'subclass' : subtype_list}) 
if args.batch:
    if args.output:
        df.to_csv(args.output)
        print('file saved')
    else:
        print('need output file path')
        parser.error('missing parameters, use -h for help')
if args.smi:
    print(df)
    
print('normal termination')       
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
