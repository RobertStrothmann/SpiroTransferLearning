from rdkit import Chem
from rdkit.Chem import AllChem

import ReadInCSV as RCSV

import math
import os

import collections
from joblib import Parallel, delayed
import random

import argparse

import time

### this script can create a number of .xyz files from a .csv REINVENT scaffold-decorator output file
### this script automatically runs in parallel when len(ScaffFolder) <  available CPUs

def GetSMIList(t,scaff,ScriptPath,SubSetSize):
    Start=time.time()
    G_CWD=os.getcwd()

    os.mkdir('{}/{}/PropPred/'.format(G_CWD,scaff))

    PathToCSV='{}/Gen_Mol.csv'.format(scaff)

        #### variable assignment
        ## PathToWorkFolder
    SplitPath=PathToCSV.split('/')[0:-1]
   
    SubSetSize=int(SubSetSize)

    WorkFolderLen=0
    for i in SplitPath:
        WorkFolderLen=WorkFolderLen+1
        WorkFolderLen=WorkFolderLen+len(i)
    
    PathToWorkFolder=PathToCSV[0:WorkFolderLen]
    
        ### Read out Different Lists of SMILES strings from .csv
    DecorationList,MoleculeList,Abundance,DecorDict=RCSV.ReadCSV(PathToCSV)
    
        #### processing smiles strings using rdkit
        ## convert string to smiles rdkit object
    Mol_list = [Chem.MolFromSmiles(smiles) for smiles in MoleculeList]
    
        #### Charge List
    ChrgList=[]
    for i in range(len(Mol_list)):
        ChrgList.append(Chem.rdmolops.GetFormalCharge(Mol_list[i]))
    
        #### Greate .xyz files in a seperate folder including weighted subset 
        #### of geometries having "SubsetSize" geometries
    
    print('step 1: after {}'.format(Start-time.time()))

        ## check for known molecules
    KnownMolFile=open('AllCalcMols.smi', 'r') 
    KnownMols=KnownMolFile.readlines()
    KnownMolFile.close()
    
    KnownCanSmiles=[]
    for i in KnownMols:
        KnownCanSmiles.append(Chem.CanonSmiles(i))
    
    
        ## created weighted mol list
    WeightedMolList=[]
    for i in range(len(MoleculeList)):
        TmpSmi=Chem.CanonSmiles(MoleculeList[i])
        
            ## kick out molecules with charged sidegroups
        if ChrgList[i] != 0:
            continue
        
            ## kick out not fully decorated molecules
        if '*' in TmpSmi:
            continue 
        
            ## kick out molecules that were already sampled
        if TmpSmi in KnownCanSmiles:
            continue
        
            ## otherwise append molecule as often as its abundance based on the sampling
        for j in range(Abundance[i]):
            WeightedMolList.append(MoleculeList[i])

    
        ## check if there are more structures available then SubSetSize variable
    WCounts=collections.Counter(WeightedMolList)
    if len(WCounts) < SubSetSize:
        SubSetSize=len(WCounts)

    print('step 2: after {}'.format(Start-time.time()))

    SubSetOfGeoms=[]
        ## create a subset of molecules to calculate
    while len(SubSetOfGeoms) < SubSetSize:
        if (len(SubSetOfGeoms)%144)==0:
            print(len(SubSetOfGeoms))

        ## discard smiles with * in it (equals only one decoration)
        TestGeomInx=random.randint(0,len(WeightedMolList)-1)
        TestSmile=WeightedMolList[TestGeomInx]
        
            ## if SMI already in subset of new geoms --> don't append
        if TestSmile in SubSetOfGeoms:
            continue

        if '*' in TestSmile:
            continue
        
            ## check if molecule is optimizable
        tmp_mol = Chem.MolFromSmiles(TestSmile,sanitize=False)
        if tmp_mol is None:
            continue

        InputMolFile=Chem.MolFromSmiles(TestSmile)
        try:
            TmpMolFile=Chem.AddHs(InputMolFile)
        except:
            print("Something went wrong")
            continue
        AllChem.AssignAtomChiralTagsFromStructure(TmpMolFile,confId=-1,replaceExistingTags=True)
        AllChem.EmbedMolecule(TmpMolFile)
        try:
            AllChem.MMFFOptimizeMolecule(TmpMolFile)
        except:
            print("Something went wrong")
        else:
            print("Nothing went wrong")
            SubSetOfGeoms.append(TestSmile)
            WeightedMolList[TestGeomInx]='*'
    
    JobList=[]
    c=0
    for mol in SubSetOfGeoms:
        JobList.append("python3 {}/SMI2XYZ.py -s {} -g {} -m \"{}\" ".format(ScriptPath,scaff,str(c).zfill(3),mol))
        c=c+1

    tmp_file=open('{}/JobList_CSV.txt'.format(scaff),'w+')
    for j in JobList:
        tmp_file.write(j)
        tmp_file.write('\n')
    tmp_file.close()

    return()

def RunParallel(t,ScriptPath,WKeys,WVals):
    ### build dict
    Weights={}
    for i in range(len(WKeys)):
        Weights[WKeys[i]]=int(WVals[i])
   
    print(Weights)

    tmp_dir=os.listdir()
    ScaffFolder=[]
    for i in tmp_dir:
        if 'Scaff_' in i:
            ScaffFolder.append(i)

    Parallel(n_jobs=len(ScaffFolder))(delayed(GetSMIList)(t,i,ScriptPath,Weights[i]) for i in ScaffFolder)

    return()

###parsing arguments
parser = argparse.ArgumentParser(description='Convert a .csv file to a list to a number of molecules.')

parser.add_argument("-t", dest="TL", required=True,
    help="TL It", metavar="TL It")

parser.add_argument("-s", dest="script", required=True,
    help="scriptpath", metavar="scriptpath")

parser.add_argument("-wk", dest="Wkeys", required=True, nargs='+', type=str,
    help="WeightsMatrixKeys", metavar="WeightsMatrixKeys")

parser.add_argument("-wv", dest="Wvals", required=True,nargs='+', type=str,
    help="WeightsMatrixVals", metavar="WeightsMatrixVals")

args = parser.parse_args()

TLIt=int(args.TL)
Scr=str(args.script)
WKeys=list(args.Wkeys)
WVals=list(args.Wvals)

## program run
RunParallel(TLIt,
          Scr,
          WKeys,
          WVals)
