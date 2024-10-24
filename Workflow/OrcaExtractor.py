import numpy as np


### reads out electronic states from Orca TDDFT
### calculation in eV, wl and F_osc
### and writes them into a file "states.txt"

def ExtractFromOrca(PathToOut,NumberOfEx=10):
    ## Define Workdir
    SplitPath=PathToOut.split('/')[0:-1]

    WorkFolderLen=0
    for i in SplitPath:
        WorkFolderLen=WorkFolderLen+1
        WorkFolderLen=WorkFolderLen+len(i)

    PathToWorkFolder=PathToOut[0:WorkFolderLen]



    ## extract states
    Schalter=0
    Zaehler=0
    State=1
    NumberOfEx=NumberOfEx+1
    ExStates=np.zeros((3,NumberOfEx-1))
    
    
    file=open('{}'.format(PathToOut),'r',errors='replace')
    for line in file:
        if line.startswith('         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS'):
            Schalter=1
            
        if Zaehler==5:
            Schalter=0
            ExStates[0,State-1]=(float(line.split()[1]))
            ExStates[1,State-1]=(float(line.split()[2]))
            ExStates[2,State-1]=(float(line.split()[3]))
            State=State+1
            
        if State==NumberOfEx:
            Zaehler=0
        
        if Schalter==1:
            Zaehler=Zaehler+1
    
    file.close()  

    OutputFile=open('{}/States.txt'.format(PathToWorkFolder),'w+')
    for line in range(NumberOfEx-1):
        OutputFile.write(str(line))
        OutputFile.write('\t')
        OutputFile.write(str(ExStates[1,line]))
        OutputFile.write('\t')
        OutputFile.write(str(ExStates[2,line]))
        OutputFile.write('\n')

    OutputFile.close()
    return()
