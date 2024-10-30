### global python packages
import os
import time
import shutil
import sys
import numpy as np
from joblib import Parallel, delayed

### make scripts from subscript folder importable
PathToScripts='/path/to/ScriptFolder'   ### change this path to a  path with all scripts
sys.path.append(PathToScripts)

### local python scripts
import HPC_Interaction
import FileManager

### part0
import Training_Part

### part1
import Sampling_Part

### part2
import RunCSV2XYZ

### part4
import GenSlurmScripts

### part5
import ParetoAnalysis

## global
PathToScaffDecModel='/path/to/reinvent-scaffold-decorator-master'    ### change this path to a path with the installed REINVT model
sys.path.append(PathToScaffDecModel)

## Training
epochs=50
batchsize_training=400

### here one has to specify a job file according to "JobFileTemplate"
JobFile=open('JobFile.txt','r')
counter=0
for line in JobFile:
    if counter==0:
        ScaffRunName=line.strip()
    if counter==1:
        Scaff_Folder=line.strip()
    counter=counter+1
JobFile.close()

## Sampling
NumOfGPUs=8

## single 256 / double 36 /triple 12
NumRand=36
NumAtt=36

## CSV2XYZ
SubSetSize=144
TLIterations=2
ThMethod='TDDFT'

## init file with all hit molecules
FileManager.InitGlobalFiles()

for t in range(TLIterations):
    if step==0:
        ## initialize Transfer learning folder
        FileManager.InitTLFolder(t)
        os.mkdir('TL_{}/TJobFiles/'.format(t))
    
        #####################################################################################################################
        ##### 1.) Training                      ---> for argparse: conditions file, baseline model, path to reinvent (global)
        start=time.time()
        Training_Part.PerformTraining(t,
                                   PathToScaffDecModel,
                                   PathToScripts,
                                  'conditions.json.example',
                                  'Baseline.model',         ### this model has to be downloaded or another pretrained model can be used
                                  epochs,
                                  batchsize_training)
    
        stop=time.time()
        FileManager.UpdateTimings('Training',round(stop - start,5))
        #####################################################################################################################
        step=step+1
    
    if step==1:
        #####################################################################################################################
        ##### 2.) Sampling Scaffolds            ---> for argparse: model for sampling, modelpath (global), 
        ##                                                              #rand. steps (global), #samples strucs (global)
        ##                                                              #GPU cores (global), list of scaffolds (global)
        start=time.time()
        ScaffFolder=Sampling_Part.PerformSampling(t,
                                              PathToScaffDecModel,
                                              PathToScripts,
                                              'Baseline.model',  ### this model has to be downloaded or another pretrained model can be used
                                              NumOfGPUs,
                                              Scaff_Folder,
                                              NumRand,
                                              NumAtt) 
      
        stop=time.time()
        FileManager.UpdateTimings('Sampling',round(stop - start,5))
        step=step+1
        #####################################################################################################################

    if step==2:
        #####################################################################################################################
        #### 3.) Create Photoswitch Subset      ---> for argparse: nothing
        ## init weights
        if t==0: 
            Weights={}
            for s in ScaffFolder:
                Weights[s]=SubSetSize
        ##################################### CPU slurm script #############################################################
        ##################################### ############### #############################################################
        start=time.time()

        ## copy all already calculated to the iteration folder
        shutil.copy('AllCalcMols.smi','TL_{}/AllCalcMols.smi'.format(t)) 

        CWD=os.getcwd()
        TL_Path='{}/TL_{}'.format(CWD,t)
        os.chdir('TL_{}'.format(t))
        ## generate Mol2XYZ input files including gnu parallel
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subCSV2SMI.sl')


        ### convert weight dict to keys and vals
        W_keys=list(Weights.keys())
        W_vals=list(Weights.values())

        W_keys_Input=''
        W_vals_Input=''
        for y in range(len(W_keys)):
            if y==len(W_keys)-1:
                W_keys_Input=W_keys_Input+str(W_keys[y])
                W_vals_Input=W_vals_Input+str(W_vals[y])
            else:
                W_keys_Input=W_keys_Input+'{} '.format(str(W_keys[y]))
                W_vals_Input=W_vals_Input+'{} '.format(str(W_vals[y]))

        tmp_file=open('subCSV2SMI.sl','a')
        tmp_file.write('python3 {}/CSV2XYZ.py -t {} -s {} -wk {} -wv {}'.format(PathToScripts,t,PathToScripts,W_keys_Input,W_vals_Input))
        tmp_file.close()

        os.system('sbatch subCSV2SMI.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'CSV2SMI')

        os.system('rm subCSV2SMI.sl')

        ### creat slurm for smi 2 xyz
        RunCSV2XYZ.CreateSubFiles(t,PathToScripts,TL_Path)

        ### run CSV2XYZ slurm files
        tmp_dir=os.listdir()
        slurm_files=[]
        cmd_files=[]
        task_files=[]
        for tmp_s_file in tmp_dir: 
            if 'subCSV2XYZ_' in tmp_s_file:
                slurm_files.append(tmp_s_file)
            if 'cmd_' in tmp_s_file:
                cmd_files.append(tmp_s_file)
        
        for slfile in slurm_files:
            os.system('sbatch {}'.format(slfile))
        
        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'xTB-TDDFT')
    
        tmp_dir=os.listdir()
        for tmp_s_file in tmp_dir:
            if 'task_' in tmp_s_file:
                task_files.append(tmp_s_file)
    
        for slfile in slurm_files:
            os.system('rm {}'.format(slfile))

        for cmfile in cmd_files:
            os.system('rm {}'.format(cmfile))

        for tfiles in task_files:
            os.system('rm {}'.format(tfiles))

        os.chdir(CWD)

        stop=time.time()
        FileManager.UpdateTimings('CreateXYZ',round(stop - start,5))
        step=step+1
        #####################################################################################################################


    if step==3:
        #####################################################################################################################
        ### 4) xTB conformer optimziation ---> for argparse: nothing
        start=time.time()
        ## xTB
        CWD=os.getcwd()

        ### gen slurm scripts and move failed runs
        CWD=os.getcwd()
        TL_Path='{}/TL_{}'.format(CWD,t)

        ### check for failed xTB runs
        os.chdir('TL_{}'.format(t))
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subGenSlurm.sl')

        tmp_file=open('subGenSlurm.sl','a')
        tmp_file.write('python3 {}/MoveFailedRDKit.py -p {} -s {}'.format(PathToScripts,TL_Path,PathToScripts))
        tmp_file.close()

        os.system('sbatch subGenSlurm.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'xTB_Slurm')

        os.system('rm subGenSlurm.sl')
        os.chdir(CWD)

        
        ### run xTB
        os.chdir('TL_{}'.format(t))
        tmp_dir=os.listdir()
        slurm_files=[]
        cmd_files=[]
        task_files=[]
        for tmp_s_file in tmp_dir:
            if 'subxTB_' in tmp_s_file:
                slurm_files.append(tmp_s_file)
            if 'cmd_' in tmp_s_file:
                cmd_files.append(tmp_s_file)

        for slfile in slurm_files:
            os.system('sbatch {}'.format(slfile))

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'xTB-TDDFT')

        tmp_dir=os.listdir()
        for tmp_s_file in tmp_dir:
            if 'task_' in tmp_s_file:
                task_files.append(tmp_s_file)

        for slfile in slurm_files:
            os.system('rm {}'.format(slfile))

        for cmfile in cmd_files:
            os.system('rm {}'.format(cmfile))

        for tfiles in task_files:
            os.system('rm {}'.format(tfiles))

        os.chdir(CWD)
        step=step+1
        stop=time.time()

        ### check for failed xTB runs
        os.chdir('TL_{}'.format(t))
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subMoveFailed.sl')

        tmp_file=open('subMoveFailed.sl','a')
        tmp_file.write('python3 {}/MoveFailedRuns.py -p {} -m {}'.format(PathToScripts,TL_Path,'xTB'))
        tmp_file.close()

        os.system('sbatch subMoveFailed.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'MoveFailed_xTB')

        os.system('rm subMoveFailed.sl')
        os.chdir(CWD)

        ### eval Conf Runs for open form 
        os.chdir('TL_{}'.format(t))
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subGenConf.sl')

        tmp_file=open('subGenConf.sl','a')
        tmp_file.write('python3 {}/GetMinConformer.py -p {}'.format(PathToScripts,TL_Path))
        tmp_file.close()

        os.system('sbatch subGenConf.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'GenConf')

        os.system('rm subGenConf.sl')
        os.chdir(CWD)

        FileManager.UpdateTimings('Crit Calc',round(stop - start,5))
        
        ### kick out wrongly converged strucs
        os.chdir('TL_{}'.format(t))
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subMoveFailedConv.sl')

        tmp_file=open('subMoveFailedConv.sl','a')
        tmp_file.write('python3 {}/CheckCoordNum.py'.format(PathToScripts))
        tmp_file.close()

        os.system('sbatch subMoveFailedConv.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'MoveUnconverged')

        os.system('rm subMoveFailedConv.sl')
        os.chdir(CWD)
        stop=time.time()
        FileManager.UpdateTimings('TDDFT Move Failed Files',round(stop - start,5))
        step=step+1
        
    #####################################################################################################################
    #### 5) Calculate spectra via TDDFT
    if step==4:
        start=time.time()
        CWD=os.getcwd()
        TL_Path='{}/TL_{}'.format(CWD,t)

        ### calc TDDFT/sTDA
        if ThMethod=='TDDFT':
            GenSlurmScripts.GenSlurmFiles_TDDFT('{}/TL_{}'.format(CWD,t),ScaffFolder,PathToScripts)

            CWD=os.getcwd()
            os.chdir('TL_{}'.format(t))
            tmp_dir=os.listdir()
            slurm_files=[]
            cmd_files=[]
            task_files=[]
            for tmp_s_file in tmp_dir:
                if 'sub_TDDFT' in tmp_s_file:
                    slurm_files.append(tmp_s_file)
                if 'cmd_' in tmp_s_file:
                    cmd_files.append(tmp_s_file)

            for slfile in slurm_files:
                os.system('sbatch {}'.format(slfile))

            time.sleep(5)
            TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

            HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'TDDFT')

            tmp_dir=os.listdir()
            for tmp_s_file in tmp_dir:
                if 'task_' in tmp_s_file:
                    task_files.append(tmp_s_file)

            for slfile in slurm_files:
                os.system('rm {}'.format(slfile))

            for cmfile in cmd_files:
                 os.system('rm {}'.format(cmfile))

            for tfiles in task_files:
                os.system('rm {}'.format(tfiles))

            os.chdir(CWD)

        ### eval failed TDDFT run
        os.chdir('TL_{}'.format(t))
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subMoveFailed.sl')

        tmp_file=open('subMoveFailed.sl','a')
        tmp_file.write('python3 {}/MoveFailedRuns.py -p {} -m {}'.format(PathToScripts,TL_Path,ThMethod))
        tmp_file.close()

        os.system('sbatch subMoveFailed.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'MoveFailed_{}'.format(ThMethod))

        os.system('rm subMoveFailed.sl')
        os.chdir(CWD)
        stop=time.time()
        FileManager.UpdateTimings('TDDFT Move Failed Files',round(stop - start,5))
        step=step+1

    #####################################################################################################################
    #### 5) Calculate Design Criteria
    if step==5:
        start=time.time()
        S_CWD=os.getcwd()
        os.chdir('TL_{}'.format(t))
        ### run crit calc slurm files
        os.chdir('TL_{}'.format(t))
        
        ### gen crit calc slurm files
        shutil.copy('{}/SlurmScript.sl'.format(PathToScripts),'subCritAnaSetup.sl')
        tmp_file=open('subCritAnaSetup.sl','a')
        tmp_file.write('python3 {}/CritAnalysis.py -s {} -t {} -m {}'.format(PathToScripts,PathToScripts,t,ThMethod))
        tmp_file.close()

        os.system('sbatch subCritAnaSetup.sl')

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        os.system('rm subCritAnaSetup.sl')

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd(),t),'CritAnaSetup')

        ### submit slurm scripts
        tmp_dir=os.listdir()
        slurm_files=[]
        for tmp_s_file in tmp_dir:
            if 'subCritAna_' in tmp_s_file:
                slurm_files.append(tmp_s_file)

        for slfile in slurm_files:
            os.system('sbatch {}'.format(slfile))

        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here

        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TJobFiles'.format(os.getcwd()),'CritAna')
        os.chdir(CWD)

        os.chdir('TL_{}'.format(t))
        tmp_dir=os.listdir()
        slurm_files=[]
        cmd_files=[]
        task_files=[]
        for tmp_s_file in tmp_dir:
            if 'subCritAna_' in tmp_s_file:
                slurm_files.append(tmp_s_file)
            if 'cmd_' in tmp_s_file:
                cmd_files.append(tmp_s_file)
            if 'task_' in tmp_s_file:
                task_files.append(tmp_s_file)

        for slfile in slurm_files:
            os.system('rm {}'.format(slfile))

        for cmfile in cmd_files:
            os.system('rm {}'.format(cmfile))

        for tfiles in task_files:
            os.system('rm {}'.format(tfiles))

        os.chdir(CWD)
        stop=time.time()
        FileManager.UpdateTimings('Crit Analysis',round(stop - start,5))
        step=step+1

    if step==6:
        start=time.time()
        S_CWD=os.getcwd()
        os.chdir('TL_{}'.format(t))
        p_start=time.time()
        Weights=ParetoAnalysis.GetParetoFront(ScaffFolder,SubSetSize)
        p_stop=time.time()
        print(round(p_stop - p_start,5))

        stop=time.time()
        os.chdir(S_CWD)
        FileManager.UpdateTimings('Pareto',round(stop - start,5))    
        step=step+1
    #####################################################################################################################    
    
    
    #####################################################################################################################
    ### 5.) Merge all scaffold .smi files to full hit .smi file for the training of the next iteration
    ###     and to over all file to collect all hits                            --> argparse: nothing

    ### check for [H+] in smiles string and delete them from FullHitMols.smi & respective Crit file
    ## find wrong smiles
    CleanedOAH=open('Cleaned_OverAH.smi','w+')
    SMIFile=open('TL_{}/Hit_Mols.txt'.format(t),'r')
    SMILines=SMIFile.readlines()
    SMIFile.close()

    WrongSMIList=[]

    for i in range(len(SMILines)):
        if '[H+]' in SMILines[i]:
            WrongSMIList.append(i)

    for i in range(len(SMILines)):
        if i in WrongSMIList:
            pass
        else:
            CleanedOAH.write(SMILines[i])

    CleanedOAH.close()

    ## write all calc mols in AllCalcMols.smi
    AllMols=open('AllCalcMols.smi','a')
    SMIFile=open('TL_{}/AllMols.smi'.format(t),'r')
    SMILines=SMIFile.readlines()
    SMIFile.close()
    
    for i in range(len(SMILines)):
        AllMols.write(SMILines[i])
    
    AllMols.close()
    step=0
    #####################################################################################################################

## post folder management
# make folder
os.mkdir('{}/'.format(ScaffRunName))

# move all TL run folders
os.system('mv TL_* {}/'.format(ScaffRunName))

# move generated Molecules
os.system('mv Cleaned_OverAH.smi {}/'.format(ScaffRunName))
os.system('mv AllCalcMols.smi {}/'.format(ScaffRunName))

# move log files
os.system('mv Timings.txt {}/'.format(ScaffRunName))
os.system('mv nohup.out {}/'.format(ScaffRunName))
os.system('mv JobFile.txt {}/'.format(ScaffRunName))
