import os
import shutil
import math

### this script creates all relevant slurm (or any other type of HPC submission) files which are calling
### created jobs in JobList_CSV.txt

def CreateSubFiles(t,ScriptPath,TL_Path):

    G_CWD=os.getcwd()

    Scaff_Folder=[]
    TL_Dir=os.listdir()
    for folder in TL_Dir:
        if 'Scaff_' in folder:
            Scaff_Folder.append(folder)

    ### get all jobs
    JobList=[]
    for scaff in Scaff_Folder:
        tmp_file=open('{}/JobList_CSV.txt'.format(scaff),'r')
        tmp_lines=tmp_file.readlines()
        tmp_file.close

        for l in tmp_lines:
            JobList.append(l)

    NumOfJobs=math.ceil(len(JobList)/1440)
    Global_Counter=0
    for slurm_job in range(NumOfJobs):
        counter=0
        shutil.copy('{}/subCSV2XYZ.sl'.format(ScriptPath),'{}/subCSV2XYZ_{}.sl'.format(TL_Path,slurm_job))

        slurm_file=open('{}/subCSV2XYZ_{}.sl'.format(TL_Path,slurm_job),'a')
        cmd_file=open('{}/cmd_{}.lst'.format(TL_Path,slurm_job),'a')
        while counter < 1440:
            if Global_Counter==len(JobList):
                break

            ### xTB TDDFT
            cmd_file.write(JobList[Global_Counter])

            Global_Counter=Global_Counter+1
            counter=counter+1

        ## only use if you want to use a slurm script to call gnu parallel
        slurm_file=open('{}/subCSV2XYZ_{}.sl'.format(TL_Path,slurm_job),'a')
        slurm_file.write('\n\n')
        slurm_file.write('parallel --delay 0.2 --joblog task_{}.log --progress -j 72 < cmd_{}.lst'.format(slurm_job,slurm_job))
        slurm_file.close()
        cmd_file.close()

    return()
