import os
import numpy as np


## helper function to initialize all necessary global files
def InitGlobalFiles():
    CleanedOAH=open('Cleaned_OverAH.smi','a')
    CleanedOAH.close()
    Tmp_File=open('AllCalcMols.smi','w+')
    Tmp_File.close()
    Timings=open('Timings.txt','a')
    Timings.close()
    
    
    
## helper function to update the timings file with a specific task and respective time
def UpdateTimings(Task,Time):
    Timings=open('Timings.txt','a')
    Timings.write('{} : {} \n'.format(Task,Time))
    Timings.close()
    
## function to create a new folder for a TL Iteration (TIt)
def InitTLFolder(TIt):
    os.mkdir('TL_{}'.format(TIt))
    Timings=open('Timings.txt','a')
    Timings.write('TL _{}:\n'.format(TIt))
    Timings.close()
