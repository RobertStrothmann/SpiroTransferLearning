from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

import ase
from ase.io import read
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
from ase.neighborlist import natural_cutoffs


### RdKit and ase script to get the coordination number
### for a given open SP molecule for specific atoms

######### ref mol 
RefSmi='C[N+]1=C(C=Cc2ccccc2[O-])C(C)(C)c2ccccc21'
RefMol=Chem.MolFromSmiles(RefSmi)

def CoordNum(PathToGeom): 
    ### conformer
    test_file=open('{}/Molecule.smi'.format(PathToGeom))
    smi=test_file.readlines()[0]
    test_file.close()

    ### load data
    TestSMI=smi
    TestMol=Chem.MolFromSmiles(TestSMI)
        
    TestXYZ='{}/xtbopt.xyz'.format(PathToGeom)  

    SMatch=TestMol.GetSubstructMatch(RefMol)

    temp_mol=read(TestXYZ)
    cut = natural_cutoffs(temp_mol,mult=1.5)
    NL=neighbor_list('i',temp_mol,cut)
    coord = np.bincount(NL)
          
    return(coord[SMatch[4]])

def CoordNumTwo(PathToGeom):
    ### conformer
    test_file=open('{}/Molecule.smi'.format(PathToGeom))
    smi=test_file.readlines()[0]
    test_file.close()

    ### load data
    TestSMI=smi
    TestMol=Chem.MolFromSmiles(TestSMI)

    TestXYZ='{}/xtbopt.xyz'.format(PathToGeom)

    SMatch=TestMol.GetSubstructMatch(RefMol)

    temp_mol=read(TestXYZ)
    cut = natural_cutoffs(temp_mol,mult=1.5)
    NL=neighbor_list('i',temp_mol,cut)
    coord = np.bincount(NL)

    return(coord[SMatch[3]])

def CoordNumOx(PathToGeom):
    ### conformer
    test_file=open('{}/Molecule.smi'.format(PathToGeom))
    smi=test_file.readlines()[0]
    test_file.close()

    ### load data
    TestSMI=smi
    TestMol=Chem.MolFromSmiles(TestSMI)

    TestXYZ='{}/xtbopt.xyz'.format(PathToGeom)

    SMatch=TestMol.GetSubstructMatch(RefMol)

    temp_mol=read(TestXYZ)
    cut = natural_cutoffs(temp_mol,mult=1.2)
    NL=neighbor_list('i',temp_mol,cut)
    coord = np.bincount(NL)

    return(coord[SMatch[11]])

def GetRing(PathToGeom):
    ### conformer
    test_file=open('{}/Molecule.smi'.format(PathToGeom))
    smi=test_file.readlines()[0]
    test_file.close()

    ### load data
    TestSMI=smi
    TestMol=Chem.MolFromSmiles(TestSMI)

    TestXYZ='{}/xtbopt.xyz'.format(PathToGeom)

    SMatch=TestMol.GetSubstructMatch(RefMol)

    TAtom=TestMol.GetAtomWithIdx(SMatch[4])
    return(TAtom.IsInRing())
