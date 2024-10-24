import shutil
import os
import argparse
import math
import sys
import time

### Script that writes a list of jobs which analyze the results of
### the TDDFT calculations according to the two properties of interes
### namely Addressability and Thermostability

def GenCritFiles(t,ScriptPath,Method='TDDFT'):
    Scaff_Folder=[]

    TL_Path=os.getcwd()

    TL_Dir=os.listdir()
    for folder in TL_Dir:
        if 'Scaff_' in folder:
            Scaff_Folder.append(folder)

    JobList=[]
    ### write job files for crit analysis
    for scaff in Scaff_Folder:
        os.chdir('{}/PropPred/'.format(scaff))
        Last_DirList=os.listdir('.')
        OutList_open=[]
        for i in Last_DirList:
            if 'O_Geom' in i:
                OutList_open.append(i)


        for geom in OutList_open:
            JobList.append("python3 {}/OpenClosedAnalysis.py -g {}/{}/PropPred/{} -m {}".format(ScriptPath,TL_Path,scaff,geom,Method))
            JobList.append("python3 {}/ThermoStabilityAnalysis.py -g {}/{}/PropPred/{}".format(ScriptPath,TL_Path,scaff,geom))

        os.chdir(TL_Path)

    NumOfJobs=math.ceil(len(JobList)/300)
    Global_Counter=0
    for slurm_job in range(NumOfJobs):
        counter=0
        shutil.copy('{}/subCritAna.sl'.format(ScriptPath),'{}/subCritAna_{}.sl'.format(TL_Path,slurm_job))

        slurm_file=open('{}/subCritAna_{}.sl'.format(TL_Path,slurm_job),'a')
        cmd_file=open('{}/cmd_{}.lst'.format(TL_Path,slurm_job),'a')
        while counter < 300:
            if Global_Counter==len(JobList):
                break

            ### xTB TDDFT
            cmd_file.write(JobList[Global_Counter])
            cmd_file.write('\n')

            Global_Counter=Global_Counter+1
            counter=counter+1

        slurm_file.close()

        ## only use if you want to use a slurm script to call gnu parallel
        slurm_file=open('{}/subCritAna_{}.sl'.format(TL_Path,slurm_job),'a')
        slurm_file.write('\n\n')
        slurm_file.write('parallel --delay 0.2 --joblog task_{}.log --progress -j 72 < cmd_{}.lst'.format(slurm_job,slurm_job))
        slurm_file.close()
        cmd_file.close()

    return()

# parsing arguments
parser = argparse.ArgumentParser(description='Create submission files for crit analysis jobs')

parser.add_argument("-s", dest="ScriptPath", required=True,
    help="Path to script folder", metavar="Path to script folder")

parser.add_argument("-t", dest="TLIt", required=True,
    help="TL iteration", metavar="TL iteration")

parser.add_argument("-m", dest="Method", required=True,
    help="QC Method", metavar="QC Method")

args = parser.parse_args()

Scripts=str(args.ScriptPath)
t=int(args.TLIt)
THMethod=str(args.Method)

## program run
GenCritFiles(t,Scripts,THMethod)
