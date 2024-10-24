from rdkit import Chem
from rdkit.Chem import AllChem

import ReadInCSV as RCSV

import os
import shutil

import collections

import random

import argparse

### script which creates a closed and open .xyz file
### for a given SP SMILES string

def CreateXYZ(SMI,Scaff,GeomID):
   
    CWD=os.getcwd()

    SubSetOfGeoms=[SMI]
    ### SMARTS SpiroOpening
    Open_SubSetGeoms=[]
    
    ### reaction
    rxn = AllChem.ReactionFromSmarts('[N:0]-[#6;D4:1]([#6:2])-[#8:3]>>([N+:0]=[#6;D4:1]([#6:2]).[#8-:3])')
    
    for i in SubSetOfGeoms:
        ### input smile
        smi1 = i
        mol1 = Chem.MolFromSmiles(smi1,sanitize=True)
        
        ### perform cleavage
        Open_SubSetGeoms.append([Chem.MolToSmiles(x) for x in rxn.RunReactants((mol1,))[0]][0])
    


    # closed
    ## make geom folder
    os.mkdir('{}/{}/PropPred/C_Geom_{}/'.format(CWD,Scaff,GeomID))
        
    ## write out smiles string
    SMIOutput=open('{}/{}/PropPred/C_Geom_{}/Molecule.smi'.format(CWD,Scaff,GeomID),'w+')
    SMIOutput.write(SubSetOfGeoms[0])
    SMIOutput.close()
   
    ## create molecule xyz coord by MMFF embedding
    InputMolFile=Chem.MolFromSmiles(SubSetOfGeoms[0])

    TmpMolFile=Chem.AddHs(InputMolFile)
    AllChem.AssignAtomChiralTagsFromStructure(TmpMolFile,confId=-1,replaceExistingTags=True)
    AllChem.EmbedMultipleConfs(TmpMolFile,numConfs=10)
    AllChem.MMFFOptimizeMolecule(TmpMolFile)
    for i in range(10):
        os.mkdir('{}/{}/PropPred/C_Geom_{}/Conf_{}'.format(CWD,Scaff,GeomID,i))
        Chem.MolToXYZFile(TmpMolFile,'{}/{}/PropPred/C_Geom_{}/Conf_{}/Input.xyz'.format(CWD,Scaff,GeomID,i),confId=i)
        SMIOutput=open('{}/{}/PropPred/C_Geom_{}/Conf_{}/Molecule.smi'.format(CWD,Scaff,GeomID,i),'w+')
        SMIOutput.write(SubSetOfGeoms[0])
        SMIOutput.close()
    
    # open
    ## make geom folder
    os.mkdir('{}/{}/PropPred/O_Geom_{}/'.format(CWD,Scaff,GeomID))
        
    ## write out smiles string
    SMIOutput=open('{}/{}/PropPred/O_Geom_{}/Molecule.smi'.format(CWD,Scaff,GeomID),'w+')
    SMIOutput.write(Open_SubSetGeoms[0])
    SMIOutput.close()
        
    ## create molecule xyz coord by MMFF embedding
    InputMolFile=Chem.MolFromSmiles(Open_SubSetGeoms[0])
    
    TmpMolFile=Chem.AddHs(InputMolFile)
    AllChem.AssignAtomChiralTagsFromStructure(TmpMolFile,confId=-1,replaceExistingTags=True)
    AllChem.EmbedMultipleConfs(TmpMolFile,numConfs=10)
    AllChem.MMFFOptimizeMolecule(TmpMolFile)
    for i in range(10):
        os.mkdir('{}/{}/PropPred/O_Geom_{}/Conf_{}'.format(CWD,Scaff,GeomID,i))
        Chem.MolToXYZFile(TmpMolFile,'{}/{}/PropPred/O_Geom_{}/Conf_{}/Input.xyz'.format(CWD,Scaff,GeomID,i),confId=i) 
        SMIOutput=open('{}/{}/PropPred/O_Geom_{}/Conf_{}/Molecule.smi'.format(CWD,Scaff,GeomID,i),'w+')
        SMIOutput.write(Open_SubSetGeoms[0])
        SMIOutput.close()

    return()
# parsing arguments 
parser = argparse.ArgumentParser(description='Convert a .csv file to a list to a number of molecules.')

parser.add_argument("-s", dest="scaff", required=True,
    help="scaffold path", metavar="scaffold path")

parser.add_argument("-m", dest="mol", required=True,
    help="Molecule as SMI", metavar="Molecule as SMI")

parser.add_argument("-g", dest="geom", required=True,
    help="Geom ID", metavar="Geom ID")

args = parser.parse_args()

SMI_Mol=str(args.mol)
Scaff=str(args.scaff)
GID=str(args.geom)

## program run
CreateXYZ(SMI_Mol,
          Scaff,
          GID)
