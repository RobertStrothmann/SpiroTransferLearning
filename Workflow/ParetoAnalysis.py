import os
import collections
import numpy as np

import argparse

import random

### script to read in all addressability and thermostability values
### for one transfer learning iteration and apply non-dominated sorting
### on it

### writes the hits in a separat file and updates the weights list

def GetParetoFront(ScaffFolder,SubSetSize,NDSLayers=3):
    Address=[]
    Thermo=[]

    Address_old=[]
    Thermo_old=[]

    Index=[]
    Scaff_ID=[]
    SMI=[]
   
    CWD=os.getcwd()

    for s in ScaffFolder:
        tmp_cwd=os.getcwd()
        os.chdir('{}/{}/PropPred/'.format(CWD,s))
        tmp_folders=os.listdir()
        
        G_IDs=[]
        for f in tmp_folders:
            if 'O_Geom' in f:
                G_IDs.append(f)

        for g in G_IDs:
            try:
                tmp_address=open('{}/{}/PropPred/{}/AddressCrit.txt'.format(CWD,s,g),'r')
            except:
                break
            
            try:
                tmp_thermo=open('{}/{}/PropPred/{}/ThermostabCrit.txt'.format(CWD,s,g),'r')
            except:
                break

            try:
                tmp_smi=open('{}/{}/PropPred/{}/Molecule.smi'.format(CWD,s,g.replace('O','C')),'r')
            except:
                break
            
            Address.append(float(tmp_address.readlines()[0]))
            Thermo.append(float(tmp_thermo.readlines()[0]))
            Index.append(s+'/{}'.format(g))
            SMI.append(tmp_smi.readlines()[0])
            Scaff_ID.append(s)

            Address_old.append(Address[-1])
            Thermo_old.append(Thermo[-1])

            tmp_address.close()
            tmp_thermo.close()
            tmp_smi.close()


        os.chdir(tmp_cwd)

    ### write out all mols
    tmp_smi=open('{}/AllMols.smi'.format(CWD),'w+')
    for mol in SMI:
        tmp_smi.write(mol)
        tmp_smi.write('\n')
    tmp_smi.close()

    
    ### define pareto points
    ParetoPoints=[]
    
    HitSMIs=[]
    HitLinks=[]
    HitAddress=[]
    HitThermo=[]
    HitScaffID=[]
    
    for n in range(NDSLayers):
        New_ParetoPoints=[]
        for i in range(len(Address)):
            Opt_Condition=True
            for j in range(len(Thermo)):
                if i==j:
                    continue
                ## check for pareto points
                if Address[i] < Address[j] and Thermo[i] < Thermo[j]:
                    Opt_Condition=False
                if Opt_Condition==False:
                    break
            if Opt_Condition==True:
                New_ParetoPoints.append(i)
                
        for p in New_ParetoPoints:
            HitSMIs.append(SMI[p])
            HitLinks.append(Index[p])
            
            HitAddress.append(Address[p])
            Address[p]=0
            
            HitThermo.append(Thermo[p])
            Thermo[p]=0
            ParetoPoints.append(p)
        
            HitScaffID.append(Scaff_ID[p])

    ### write files
    tmp_Addressfile=open('{}/Hit_Address.txt'.format(CWD),'w+')
    tmp_Thermofile=open('{}/Hit_Thermo.txt'.format(CWD),'w+')
    tmp_Indexfile=open('{}/Hit_Index.txt'.format(CWD),'w+')
    tmp_Molsfile=open('{}/Hit_Mols.txt'.format(CWD),'w+')

    for i in range(len(HitSMIs)):
        tmp_Addressfile.write(str(HitAddress[i]))
        tmp_Addressfile.write('\n')
        
        tmp_Thermofile.write(str(HitThermo[i]))
        tmp_Thermofile.write('\n')
        
        tmp_Indexfile.write(str(HitLinks[i]))
        tmp_Indexfile.write('\n')
        
        tmp_Molsfile.write(str(HitSMIs[i]))
        tmp_Molsfile.write('\n')
        
    tmp_Addressfile.close()
    tmp_Thermofile.close()
    tmp_Indexfile.close()
    tmp_Molsfile.close()
    

    #### define weights
    Weights={}

    FullNumStrucs=len(ScaffFolder)*SubSetSize

    ### get collection of IndexList
    HitIndexCollection=collections.Counter(HitScaffID)

    for w in ScaffFolder:
        Weights[w]=5+int(HitIndexCollection[w]/sum(list(HitIndexCollection.values()))*FullNumStrucs)

    print(Weights)
    return(Weights)


