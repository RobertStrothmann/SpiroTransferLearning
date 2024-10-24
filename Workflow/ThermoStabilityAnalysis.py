from rdkit import Chem

import numpy as np
import os
import shutil

import argparse

import sys

### script to evaluate a run folder for a given geometry
### according to the thermostability criterion based on the open form

### writes a file "ThermostabCrit.txt" with the result in that folder

def AnalyzeThermostab(GeomPath):
    ### Global variables
    RefSmi='C[N+]1=C(/C=C/C2=CC=CC=C2[O-])C(C)(C)C3=C1C=CC=C3'
    RefMol=Chem.MolFromSmiles(RefSmi)
    
    GID=GeomPath.split('/')[-1]
    TreePath=''
    for j in GeomPath.split('/')[0:-1]:
        TreePath=TreePath+j+'/'

    ### Spiro C index: 2
    ### Spiro O index: 11
    ### DB C index1: 3
    ### DB C index2: 4

    tmp_file=open('{}/Molecule.smi'.format(GeomPath),'r')
    tmp_smi=tmp_file.readlines()[0]
    mol=Chem.MolFromSmiles(tmp_smi)

    coords=np.genfromtxt('{}/xtbopt.xyz'.format(GeomPath),skip_header=2)
    coords=coords[:,1:]

    MatchList=mol.GetSubstructMatch(RefMol)

    if len(MatchList) < 1:
        tmp_file_address=open('{}/ThermostabCrit.txt'.format(GeomPath),'w+')
        tmp_file_address.write(str(0))
        tmp_file_address.close()
        shutil.move('{}/'.format(GeomPath),'{}/Failed/{}'.format(TreePath,GID))
        shutil.move('{}/'.format(GeomPath.replace('O_','C_')),'{}/Failed/{}'.format(TreePath,GID.replace('O_','C_')))
        return()

    if len(coords) < 1:
        tmp_file_address=open('{}/ThermostabCrit.txt'.format(GeomPath),'w+')
        tmp_file_address.write(str(0))
        tmp_file_address.close()
        shutil.move('{}/'.format(GeomPath),'{}/Failed/{}'.format(TreePath,GID))
        shutil.move('{}/'.format(GeomPath.replace('O_','C_')),'{}/Failed/{}'.format(TreePath,GID.replace('O_','C_')))
        return()

    ### check if it is really open - if O spiro and
    try:
        coords[MatchList[11]]
    except:
        tmp_file_address=open('{}/ThermostabCrit.txt'.format(GeomPath),'w+')
        tmp_file_address.write(str(0))
        tmp_file_address.close()
        shutil.move('{}/'.format(GeomPath),'{}/Failed/{}'.format(TreePath,GID))
        shutil.move('{}/'.format(GeomPath.replace('O_','C_')),'{}/Failed/{}'.format(TreePath,GID.replace('O_','C_')))
        return()

    if np.linalg.norm(coords[MatchList[2]]-coords[MatchList[11]]) < 2:
        tmp_file_address=open('{}/ThermostabCrit.txt'.format(GeomPath),'w+')
        tmp_file_address.write(str(0))
        tmp_file_address.close()
        shutil.move('{}/'.format(GeomPath),'{}/Failed/{}'.format(TreePath,GID))
        shutil.move('{}/'.format(GeomPath.replace('O_','C_')),'{}/Failed/{}'.format(TreePath,GID.replace('O_','C_')))

        
    else:
        TmpRes=np.linalg.norm(coords[MatchList[3]]-coords[MatchList[4]])
        tmp_file_address=open('{}/ThermostabCrit.txt'.format(GeomPath),'w+')
        tmp_file_address.write(str(TmpRes))
        tmp_file_address.close()
            
    return()

# parsing arguments 
parser = argparse.ArgumentParser(description='Calculate the thermostability value for value for a folder of geom folders')

parser.add_argument("-g", dest="GePath", required=True,
    help="Path to geom folders", metavar="Path to geom folders")

args = parser.parse_args()

InpGeom=str(args.GePath)

## program run
AnalyzeThermostab(InpGeom)

