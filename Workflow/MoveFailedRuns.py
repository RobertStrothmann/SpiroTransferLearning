import shutil
import os
import argparse
import math
import sys
import time
from joblib import Parallel, delayed

### this script helps to move failed xTB geometry optimization
### and failed orca TDDFT calculations to a "FailedRuns" directory

### in the end of the script there is a implementet paralelization
### for len(ScaffFolder) cores

def MoveFailedRuns_xTB(scaff,TL_Path):
    os.chdir('{}/{}/PropPred/'.format(TL_Path,scaff))

    GeomPath=os.getcwd()
    #### check for converged xTB
    DirList=os.listdir()
    OutList_closed=[]
    for i in DirList:
        if 'C_Geom' in i:
            OutList_closed.append(i)

    OutList_open=[]
    for i in DirList:
        if 'O_Geom' in i:
            OutList_open.append(i)

    FailedOut_open=[]
    FailedOut_closed=[]

    for i in range(len(OutList_closed)):
        print('___')
        print(scaff)
        print(i)
        print('___')
        BreakSwitch=False
        for c in range(10):
            if BreakSwitch==True:
                break
            try:
                TmpOutFile=open('{}/{}/Conf_{}/out'.format(GeomPath,OutList_closed[i],c),errors='replace')
            except:
                FailedOut_closed.append(OutList_closed[i])
                BreakSwitch=True
                break

            TmpOutFile=open('{}/{}/Conf_{}/out'.format(GeomPath,OutList_closed[i],c),errors='replace')
            OutFileLines=TmpOutFile.readlines()
            TmpOutFile.close()

            for line in OutFileLines[-20:]:
                if '[ERROR] Program stopped' in line:
                    BreakSwitch=True
                    FailedOut_closed.append(OutList_closed[i])

    for i in range(len(OutList_open)):
        print('___')
        print(scaff)
        print(i)
        print('___')
        BreakSwitch=False
        for c in range(10):
            if BreakSwitch==True:
                break
            try:
                TmpOutFile=open('{}/{}/Conf_{}/out'.format(GeomPath,OutList_open[i],c),errors='replace')
            except:
                FailedOut_open.append(OutList_open[i])
                BreakSwitch=True
                break

            TmpOutFile=open('{}/{}/Conf_{}/out'.format(GeomPath,OutList_open[i],c),errors='replace')
            OutFileLines=TmpOutFile.readlines()
            TmpOutFile.close()

            for line in OutFileLines[-20:]:
                if '[ERROR] Program stopped' in line:
                    FailedOut_open.append(OutList_open[i])
                    BreakSwitch=True
                    break

    ## delete failed Geoms
    for i in FailedOut_open:
        shutil.move('{}/{}'.format(GeomPath,i),'{}/FailedRuns/{}'.format(GeomPath,i))
        shutil.move('{}/{}'.format(GeomPath,i.replace('O_','C_')),'{}/FailedRuns/{}'.format(GeomPath,i.replace('O_','C_')))
        if i in FailedOut_closed:
            FailedOut_closed.remove(i.replace('O_','C_'))

    for i in FailedOut_closed:
        shutil.move('{}/{}'.format(GeomPath,i),'{}/FailedRuns/{}'.format(GeomPath,i))
        shutil.move('{}/{}'.format(GeomPath,i.replace('C_','O_')),'{}/FailedRuns/{}'.format(GeomPath,i.replace('C_','O_')))

    return('{} done'.format(scaff))

def MoveFailedRuns_TDDFT(scaff,TL_Path):
    os.chdir('{}/{}/PropPred/'.format(TL_Path,scaff))

    GeomPath=os.getcwd()

    #### check for failed TDDFT
    DirList=os.listdir()
    OutList_closed=[]
    for i in DirList:
        if 'C_Geom' in i:
            OutList_closed.append(i)

    OutList_open=[]
    for i in DirList:
        if 'O_Geom' in i:
            OutList_open.append(i)

    FailedOut_open=[]
    FailedOut_closed=[]

    for i in range(len(OutList_closed)):
        try:
            TmpOutFile=open('{}/{}/TDDFT/OrcaTDDFT.out'.format(GeomPath,OutList_closed[i]),errors='replace')
        except:
            FailedOut_closed.append(OutList_closed[i])
            break

        TmpOutFile=open('{}/{}/TDDFT/OrcaTDDFT.out'.format(GeomPath,OutList_closed[i]),errors='replace')
        OutFileLines=TmpOutFile.readlines()
        TmpOutFile.close()

        tmp_lin=OutFileLines[-1]
        if 'TOTAL RUN TIME' not in tmp_lin:
            FailedOut_closed.append(OutList_closed[i])

    for i in range(len(OutList_open)):
        try:
            TmpOutFile=open('{}/{}/TDDFT/OrcaTDDFT.out'.format(GeomPath,OutList_open[i]),errors='replace')
        except:
            FailedOut_open.append(OutList_open[i])
            break

        TmpOutFile=open('{}/{}/TDDFT/OrcaTDDFT.out'.format(GeomPath,OutList_open[i]),errors='replace')
        OutFileLines=TmpOutFile.readlines()
        TmpOutFile.close()

        tmp_lin=OutFileLines[-1]
        if 'TOTAL RUN TIME' not in tmp_lin:
            FailedOut_open.append(OutList_open[i])

    ## delete failed Geoms
    for i in FailedOut_open:
        shutil.move('{}/{}'.format(GeomPath,i),'{}/Failed/{}'.format(GeomPath,i))
        shutil.move('{}/{}'.format(GeomPath,i.replace('O_','C_')),'{}/Failed/{}'.format(GeomPath,i.replace('O_','C_')))
        if i.replace('O_','C_') in FailedOut_closed:
            FailedOut_closed.remove(i.replace('O_','C_'))

    for i in FailedOut_closed:
        shutil.move('{}/{}'.format(GeomPath,i),'{}/Failed/{}'.format(GeomPath,i))
        shutil.move('{}/{}'.format(GeomPath,i.replace('C_','O_')),'{}/Failed/{}'.format(GeomPath,i.replace('C_','O_')))

    return('{} done'.format(scaff))

def MoveFailedRuns_sTDA(scaff):
    TL_Path=os.getcwd()

    os.chdir('{}/{}/PropPred/'.format(TL_Path,scaff))

    GeomPath=os.getcwd()

    #### check for failed sTDA
    DirList=os.listdir()
    OutList_closed=[]
    for i in DirList:
        if 'C_Geom' in i:
            OutList_closed.append(i)

    OutList_open=[]
    for i in DirList:
        if 'O_Geom' in i:
            OutList_open.append(i)

    FailedOut_open=[]
    FailedOut_closed=[]

    for i in range(len(OutList_closed)):
        try:
            TmpOutFile=open('{}/{}/sTDA/stda.out'.format(GeomPath,OutList_closed[i]),errors='replace')
        except:
            FailedOut_closed.append(OutList_closed[i])
            break

        TmpOutFile=open('{}/{}/sTDA/stda.out'.format(GeomPath,OutList_closed[i]),errors='replace')
        OutFileLines=TmpOutFile.readlines()
        TmpOutFile.close()

        tmp_lin=OutFileLines[-2]
        if ' sTDA done.' not in tmp_lin:
            FailedOut_closed.append(OutList_closed[i])

    for i in range(len(OutList_open)):
        try:
            TmpOutFile=open('{}/{}/sTDA/stda.out'.format(GeomPath,OutList_open[i]),errors='replace')
        except:
            FailedOut_open.append(OutList_open[i])
            break

        TmpOutFile=open('{}/{}/sTDA/stda.out'.format(GeomPath,OutList_open[i]),errors='replace')
        OutFileLines=TmpOutFile.readlines()
        TmpOutFile.close()

        tmp_lin=OutFileLines[-2]
        if ' sTDA done.' not in tmp_lin:
            FailedOut_open.append(OutList_open[i])

    ## delete failed Geoms
    for i in FailedOut_open:
        shutil.move('{}/{}'.format(GeomPath,i),'{}/Failed/{}'.format(GeomPath,i))
        shutil.move('{}/{}'.format(GeomPath,i.replace('O_','C_')),'{}/Failed/{}'.format(GeomPath,i.replace('O_','C_')))
        if i.replace('O_','C_') in FailedOut_closed:
            FailedOut_closed.remove(i.replace('O_','C_'))

    for i in FailedOut_closed:
        shutil.move('{}/{}'.format(GeomPath,i),'{}/Failed/{}'.format(GeomPath,i))
        shutil.move('{}/{}'.format(GeomPath,i.replace('C_','O_')),'{}/Failed/{}'.format(GeomPath,i.replace('C_','O_')))

    return('{} done'.format(scaff))

def SubMoveFailed(TL_Path,THMethod):
    ### get scaff folder
    ScaffFolder=[]
    tmp_dir=os.listdir()
    for i in tmp_dir:
        if 'Scaff_' in i:
            ScaffFolder.append(i)
    
    Res=[]
    for scaff in ScaffFolder:
        if THMethod=='xTB':
            Res.append(Parallel(n_jobs=len(ScaffFolder))(delayed(MoveFailedRuns_xTB)(i,TL_Path) for i in ScaffFolder))
        if THMethod=='TDDFT':
            Res.append(Parallel(n_jobs=len(ScaffFolder))(delayed(MoveFailedRuns_TDDFT)(i,TL_Path) for i in ScaffFolder))

    print(Res)

# parsing arguments
parser = argparse.ArgumentParser(description='Create submission files to move failed jobs')

parser.add_argument("-p", dest="TLPath", required=True,
    help="TLPath", metavar="TLPath")

parser.add_argument("-m", dest="Method", required=True,
    help="Method", metavar="Method")

args = parser.parse_args()

TL_Path=str(args.TLPath)
THMethod=str(args.Method)

## program run
SubMoveFailed(TL_Path,THMethod)
