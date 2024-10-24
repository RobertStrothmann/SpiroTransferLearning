import os
import numpy as np
import argparse
import shutil

from joblib import Parallel, delayed

### this script analyzes which of the 10 initially generated RdKit geometries
### is the lowest in xTB energy after xTB optimization

def GetxTBOpt(scaff,TL_Path):
    
    AllDir=os.listdir('{}/{}/PropPred/'.format(TL_Path,scaff))

    GeomFile=[]

    for file in AllDir:
        if '_Geom_' in file:
            GeomFile.append(file)


    for g in GeomFile:
        ### try if Conf folders got zipped already
        try:
            tmp_file=open('{}/{}/PropPred/{}/Conf_0/xtbopt.xyz'.format(TL_Path,scaff,g),'r')
        except:
            continue

        tmp_Ene=[]
        for c in range(10):
            tmp_file=open('{}/{}/PropPred/{}/Conf_{}/xtbopt.xyz'.format(TL_Path,scaff,g,c),'r')
            tmp_lines=tmp_file.readlines()
            tmp_file.close()

            for l in tmp_lines:
                if 'energy:' in l:
                    tmp_Ene.append(float(l.split()[1]))
                    break

        MinIdx=np.argmin(tmp_Ene)

        shutil.copy('{}/{}/PropPred/{}/Conf_{}/xtbopt.xyz'.format(TL_Path,scaff,g,MinIdx),'{}/{}/PropPred/{}/xtbopt.xyz'.format(TL_Path,scaff,g))
        shutil.copy('{}/{}/PropPred/{}/Conf_{}/out'.format(TL_Path,scaff,g,MinIdx),'{}/{}/PropPred/{}/out'.format(TL_Path,scaff,g))

        ### move and zip other confs
        os.mkdir('{}/{}/PropPred/{}/HighEConfs'.format(TL_Path,scaff,g))
        for c in range(10):
            shutil.move('{}/{}/PropPred/{}/Conf_{}/'.format(TL_Path,scaff,g,c),'{}/{}/PropPred/{}/HighEConfs/Conf_{}'.format(TL_Path,scaff,g,c))

        Bash_Comm='tar -czvf {}/{}/PropPred/{}/HighEConfs.tar.gz {}/{}/PropPred/{}/HighEConfs/'.format(TL_Path,scaff,g,TL_Path,scaff,g)
        os.system(Bash_Comm)

        Bash_Comm='rm -r {}/{}/PropPred/{}/HighEConfs/'.format(TL_Path,scaff,g)
        os.system(Bash_Comm)        

    return()


def SubGetMinConf(TL_Path):

    ### get scaff folder
    ScaffFolder=[]
    tmp_dir=os.listdir()
    for i in tmp_dir:
        if 'Scaff_' in i:
            ScaffFolder.append(i)

    for scaff in ScaffFolder:
        Parallel(n_jobs=len(ScaffFolder))(delayed(GetxTBOpt)(i,TL_Path) for i in ScaffFolder)

# parsing arguments
parser = argparse.ArgumentParser(description='Get xTB min energy conformer for open and closed')

parser.add_argument("-p", dest="TLPath", required=True,
    help="TLPath", metavar="TLPath")

args = parser.parse_args()

TL_Path=str(args.TLPath)

## program run
SubGetMinConf(TL_Path)
