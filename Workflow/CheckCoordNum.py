import numpy as np
import os
import shutil

import GetCoordNumber
from joblib import Parallel, delayed

### this script uses a RdKit routine to kick out wrongly converged spiropyran
### structures. Wrongly converged means that the open form collapsed into
### the closed form during the xTB optimization. Another way of wrongly
### converged structures are intermolecular bridges in the open form,
### which affects the hydrogen atoms at the open form (MC) rotatable double bond
### The wrong structures are moved to "FailedRuns"


def CoordCheck(s):
    CWD=os.getcwd()
    GCount=0
    os.chdir('{}/PropPred'.format(s))

    Geoms=[]
    Dir=os.listdir()
    for d in Dir:
        if 'O_Geom' in d:
            Geoms.append(d)

    for g in Geoms:
        try: 
            GetCoordNumber.CoordNum('{}/{}/PropPred/{}'.format(CWD,s,g))
        except:
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g))
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g.replace('O_','C_')))
            if g == Geoms[-1]:
                os.chdir(CWD)
                continue
            continue

        try: 
            GetCoordNumber.CoordNumTwo('{}/{}/PropPred/{}'.format(CWD,s,g))
        except:
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g))
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g.replace('O_','C_')))
            if g == Geoms[-1]:
                os.chdir(CWD)
                continue
            continue


        try: 
            GetCoordNumber.CoordNumOx('{}/{}/PropPred/{}'.format(CWD,s,g))
        except:
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g))
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g.replace('O_','C_')))
            if g == Geoms[-1]:
                os.chdir(CWD)
                continue
            continue

        if GetCoordNumber.CoordNum('{}/{}/PropPred/{}'.format(CWD,s,g)) > 3:
            GCount=GCount+1
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g))
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g.replace('O_','C_')))
            if g == Geoms[-1]:
                os.chdir(CWD)
                continue

            continue

        if GetCoordNumber.CoordNumTwo('{}/{}/PropPred/{}'.format(CWD,s,g)) > 3:
            GCount=GCount+1
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g))
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g.replace('O_','C_')))
            if g == Geoms[-1]:
                os.chdir(CWD)
                continue

            continue

        if GetCoordNumber.CoordNumOx('{}/{}/PropPred/{}'.format(CWD,s,g)) > 1:
            GCount=GCount+1
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g))
            shutil.move('{}/{}/PropPred/{}'.format(CWD,s,g.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(CWD,s,g.replace('O_','C_')))
            if g == Geoms[-1]:
                os.chdir(CWD)
                continue

            continue


    os.chdir(CWD)
    return(GCount)

ScaffFolder=[]
Dir=os.listdir()
for d in Dir:
    if 'Scaff_' in d:
        ScaffFolder.append(d)


U=Parallel(n_jobs=len(ScaffFolder))(delayed(CoordCheck)(i) for i in ScaffFolder)

print(len(U))
print(U)

