import shutil
import os
import time
import HPC_Interaction



### this script is handling the submission of scripts to query a model and
### run the sampling on a HPC system
### This function returns the names of the generated Folders including the
### sampled molecules as a .csv file

def PerformSampling(t,
                    Model_Path,
                    ScriptFolder_Path,
                    PriorModel,
                    NumOfGPU,
                    ScaffFolder,
                    NumRand,
                    NumAtt):
    
    ## variables
    CWD=os.getcwd()
    
    
    ## make a new Sampling folder for the respective transfer iteration run
    os.mkdir('TL_{}/Sampling'.format(t)) 
    
    
    ## copy slurm script and python script to the correct folder
    shutil.copy('{}/SamplingScaffold.sl'.format(ScriptFolder_Path),'{}/SamplingScaffold.sl'.format(Model_Path))
    shutil.copy('{}/RunSampling.py'.format(ScriptFolder_Path),'{}/RunSampling.py'.format(Model_Path))
    
    
    ## Write as many slurm submittion scripts as possible based on NumOfGPU
    for gpu in range(NumOfGPU):
        shutil.copy('{}/SamplingScaffold.sl'.format(Model_Path),'{}/SamplingScaffold_{}.sl'.format(Model_Path,gpu))
        SampleJobFile=open('{}/SamplingScaffold_{}.sl'.format(Model_Path,gpu),'a')
        # if t=0 the prior model is used for sampling
        if t ==0: 
            SampleJobFile.write('python3 RunSampling.py -m {}/trained_models/{} -i {} -gpu {} -num_gpu {} -scaffFolder {} -r {} -n {}'.format(
                Model_Path,PriorModel,
                Model_Path,
                NumOfGPU,
                gpu,
                ScaffFolder,
                NumRand,
                NumAtt))
    
        else:
            SampleJobFile.write('python3 RunSampling.py -m {}/TL_{}/NewModel/models/model.trained.50 -i {} -gpu {} -num_gpu {} -scaffFolder {} -r {} -n {}'.format(
                CWD,t,
                Model_Path,
                NumOfGPU,
                gpu,
                ScaffFolder,
                NumRand,
                NumAtt))
          
        SampleJobFile.close()
    
    
    ## Submitting the Job File to the Cluster
    os.chdir('{}'.format(Model_Path))
    for i in range(NumOfGPU):
        os.system('sbatch SamplingScaffold_{}.sl'.format(i))
    time.sleep(5)
    
    # move job files
    TJobNames=HPC_Interaction.AllFinished($USERNAME) ### add your HPC username here
    
    # remove not needed scripts
    for i in range(NumOfGPU):
        os.remove('SamplingScaffold_{}.sl'.format(i))
    os.remove('SamplingScaffold.sl')
    
    HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TL_{}/Sampling'.format(CWD,t),'Sampling')
    
    ##list all csv files
    CSVFiles=[]
    for i in os.listdir():
        if '.csv' in i:
            CSVFiles.append(i)
                
    #create folders for all Scaffolds
    ScaffFolder=[]
    for i in CSVFiles:
        ScaffSpecif=i.split('.')[0]
        ScaffFolder.append(ScaffSpecif) ## names based on scaffold id's
        os.mkdir('{}/TL_{}/{}'.format(CWD,t,ScaffSpecif))
        shutil.move(i,'{}/TL_{}/{}/Gen_Mol.csv'.format(CWD,t,ScaffSpecif))

    os.chdir(CWD)
    return(ScaffFolder)
