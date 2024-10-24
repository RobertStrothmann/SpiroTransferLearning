import numpy as np
import os

import OrcaExtractor as OrcEx
import PhotoswitchCriterion as PhoCrit
import shutil

import argparse

import sys

### I/O script to pass the states to the PhotoswitchCriterion script

def OpenClosedAna(GeomPath,Method='TDDFT'):
    
    if Method=='TDDFT':
        ## Extract first osc.strength and state energies (based on the init of output list)
        OrcEx.ExtractFromOrca('{}/TDDFT/OrcaTDDFT.out'.format(GeomPath.replace('O_','C_')))
        SpecA=np.genfromtxt('{}/TDDFT/States.txt'.format(GeomPath.replace('O_','C_')))

        ## Extract open osc.strength and state energies 
        OrcEx.ExtractFromOrca('{}/TDDFT/OrcaTDDFT.out'.format(GeomPath))
        SpecB=np.genfromtxt('{}/TDDFT/States.txt'.format(GeomPath))

    ### Analyse the bin spectra (open and closed form) according to the photoswitch criterion
    TmpRes=PhoCrit.PhotoswitchCrit(SpecA,SpecB)
    tmp_file_address=open('{}/AddressCrit.txt'.format(GeomPath),'w+')
    tmp_file_address.write(str(TmpRes))
    tmp_file_address.close()
        
    return()


# parsing arguments 
parser = argparse.ArgumentParser(description='Calculate the addressability value for a folder of geom folders')

parser.add_argument("-g", dest="GeomPath", required=True,
    help="Path to geom folders", metavar="Path to geom folders")

parser.add_argument("-m", dest="Method", required=True,
    help="QC Method", metavar="QC Method")

args = parser.parse_args()

InpGeom=str(args.GeomPath)
THMethod=str(args.Method)

## program run
OpenClosedAna(InpGeom,THMethod)

