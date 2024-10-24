import os
import shutil
import math


### this script generates a job list to submit to a HPC system
### in the job a new folder (TDDFT) is created and a ORCA input
### file is written

def GenSlurmFiles_TDDFT(TL_Path,ScaffFolder,ScriptPaths):
    JobList=[]

    for scaff in ScaffFolder:
        AllDir=os.listdir('{}/{}/PropPred/'.format(TL_Path,scaff))
        GeomFile=[]

        for file in AllDir:
            if '_Geom_' in file:
                GeomFile.append(file)

        for geom in GeomFile:
            Cmd_string=''

            ## TDDFT part
            Cmd_string=Cmd_string+"mkdir {}/PropPred/{}/TDDFT;".format(scaff,geom)
            Cmd_string=Cmd_string+"cp {}/OrcaTemplateTDDFT {}/PropPred/{}/TDDFT/OrcaTDDFT.inp;".format(ScriptPaths,scaff,geom)
            Cmd_string=Cmd_string+"HEADINT=$(head -n 1 {}/PropPred/{}/xtbopt.xyz);".format(scaff,geom)
            Cmd_string=Cmd_string+"tail -n $HEADINT {}/PropPred/{}/xtbopt.xyz >> {}/PropPred/{}/TDDFT/OrcaTDDFT.inp;".format(scaff,geom,scaff,geom)
            Cmd_string=Cmd_string+"echo '*' >> {}/PropPred/{}/TDDFT/OrcaTDDFT.inp;".format(scaff,geom)
            Cmd_string=Cmd_string+"cd {}/PropPred/{}/TDDFT;".format(scaff,geom)
            Cmd_string=Cmd_string+"/path/to/orca OrcaTDDFT.inp > OrcaTDDFT.out;"
            Cmd_string=Cmd_string+"cd .. ; cd .. ; cd .. ; cd .. "

            JobList.append(Cmd_string)


    NumOfJobs=math.ceil(len(JobList)/72)
    Global_Counter=0
    for slurm_job in range(NumOfJobs):
        counter=0
        shutil.copy('{}/subTDDFT_Template.sl'.format(ScriptPaths),'{}/sub_TDDFT_{}.sl'.format(TL_Path,slurm_job))

        slurm_file=open('{}/sub_TDDFT_{}.sl'.format(TL_Path,slurm_job),'a')
        cmd_file=open('{}/cmd_{}.lst'.format(TL_Path,slurm_job),'a')
        while counter < 72:
            if Global_Counter==len(JobList):
                break

            ### xTB TDDFT
            cmd_file.write(JobList[Global_Counter])
            cmd_file.write('\n')

            Global_Counter=Global_Counter+1
            counter=counter+1

        ## only use if you want to use a slurm script to call gnu parallel
        slurm_file=open('{}/sub_TDDFT_{}.sl'.format(TL_Path,slurm_job),'a')
        slurm_file.write('\n\n')
        slurm_file.write('parallel --delay 0.2 --joblog task.log --progress -j 72 < cmd_{}.lst'.format(slurm_job))
        slurm_file.close()
        cmd_file.close()
    return()
