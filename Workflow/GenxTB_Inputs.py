import os
import shutil
import math
from joblib import Parallel, delayed
import argparse

### this script evaluates the created .xyz files based on RdKit and
### purges wrong/not existing folders into a "FailedRuns" directory

### If the .xyz exists there will be a xtb input command written
### into a list of commands

def EvalScaff(scaff,TL_Path):
    JobList=[]
    AllDir=os.listdir('{}/{}/PropPred/'.format(TL_Path,scaff))
    GeomFile=[]
        
    os.mkdir('{}/{}/PropPred/FailedRuns'.format(TL_Path,scaff))
        
    for file in AllDir:
        if '_Geom_' in file:
            GeomFile.append(file)
        
    for geom in GeomFile:
        if 'C_Geom' in geom:
            for c in range(10):
                Cmd_string=''
                try:
                    with open('{}/{}/PropPred/{}/Conf_{}/Input.xyz'.format(TL_Path,scaff,geom,c)) as f:
                        HeadInt = int(f.readlines()[0])
                except:
                    try:
                        shutil.move('{}/{}/PropPred/{}'.format(TL_Path,scaff,geom),'{}/{}/PropPred/FailedRuns/{}'.format(TL_Path,scaff,geom))
                    except:
                        pass
                    try:
                        shutil.move('{}/{}/PropPred/{}'.format(TL_Path,scaff,geom.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(TL_Path,scaff,geom.replace('O_','C_')))
                    except:
                        pass
                    continue

                
                Cmd_string=Cmd_string+"cd {}/PropPred/{}/Conf_{};".format(scaff,geom,c)
                Cmd_string=Cmd_string+"/path/to/xtb Input.xyz --alpb acetonitrile --opt  >> out;"
                Cmd_string=Cmd_string+" cd .. ; cd .. ; cd .. ;cd ..;"

                JobList.append(Cmd_string)
                
        if 'O_Geom' in geom:
            for c in range(10):
                Cmd_string=''
                try:
                    with open('{}/{}/PropPred/{}/Conf_{}/Input.xyz'.format(TL_Path,scaff,geom,c)) as f:
                        HeadInt = int(f.readlines()[0])
                except:
                    try:
                        shutil.move('{}/{}/PropPred/{}'.format(TL_Path,scaff,geom),'{}/{}/PropPred/FailedRuns/{}'.format(TL_Path,scaff,geom))
                    except:
                        pass
                    try:
                        shutil.move('{}/{}/PropPred/{}'.format(TL_Path,scaff,geom.replace('O_','C_')),'{}/{}/PropPred/FailedRuns/{}'.format(TL_Path,scaff,geom.replace('O_','C_')))
                    except:
                        pass
                    continue


                    
                Cmd_string=Cmd_string+"cd {}/PropPred/{}/Conf_{}/;".format(scaff,geom,c)
                Cmd_string=Cmd_string+"/path/to/xtb Input.xyz --alpb acetonitrile --opt  >> out;"
                Cmd_string=Cmd_string+" cd .. ;cd ..; cd .. ; cd .. ;"
                
                JobList.append(Cmd_string)
    return(JobList)


def GenSlurmFiles_xTB(TL_Path,ScriptPaths):
    ### get scaff folder
    ScaffFolder=[]
    tmp_dir=os.listdir()
    for i in tmp_dir:
        if 'Scaff_' in i:
            ScaffFolder.append(i)

    JobList=[]
    
    JobRes=Parallel(n_jobs=len(ScaffFolder))(delayed(EvalScaff)(i,TL_Path) for i in ScaffFolder)

    for l in JobRes:
        JobList=JobList+l

    NumOfJobs=math.ceil(len(JobList)/1440)
    Global_Counter=0
    for slurm_job in range(NumOfJobs):
        counter=0
        shutil.copy('{}/subxTB.sl'.format(ScriptPaths),'{}/subxTB_{}.sl'.format(TL_Path,slurm_job))

        slurm_file=open('{}/subxTB_{}.sl'.format(TL_Path,slurm_job),'a')
        cmd_file=open('{}/cmd_{}.lst'.format(TL_Path,slurm_job),'a')
        while counter < 1440:
            if Global_Counter==len(JobList):
                break
            
            ### xTB TDDFT
            cmd_file.write(JobList[Global_Counter])
            cmd_file.write('\n')

            Global_Counter=Global_Counter+1
            counter=counter+1

        ## only use if you want to use a slurm script to call gnu parallel
        slurm_file=open('{}/subxTB_{}.sl'.format(TL_Path,slurm_job),'a')
        slurm_file.write('\n\n')
        slurm_file.write('parallel --delay 0.2 --joblog task_{}.log --progress -j 72 < cmd_{}.lst'.format(slurm_job,slurm_job))
        slurm_file.close()
        cmd_file.close()
    return()


# parsing arguments
parser = argparse.ArgumentParser(description='Create submission files to gen slurm scripts')

parser.add_argument("-p", dest="TLPath", required=True,
    help="TLPath", metavar="TLPath")

parser.add_argument("-s", dest="Scripts", required=True,
    help="Scripts", metavar="Scripts")

args = parser.parse_args()

TL_Path=str(args.TLPath)
ScriptsPath=str(args.Scripts)

## program run
GenSlurmFiles_xTB(TL_Path,ScriptsPath)
