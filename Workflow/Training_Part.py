import os
import shutil
import time


import HPC_Interaction


### this script uses code from REINVENT to train a model based on a set of structures
### for transfer learning iteration 0 this part of the code is not applied

def PerformTraining(t,
                    Model_Path,
                    ScriptFolder_Path,
                    Cond_File,
                    PriorModel,
                    epochs=50,
                    batchsize=1600):
    
    if t==0:
        ## for the 0 iteration no training needs to be performed and a prior model is used
        pass
    else:
        ## information on current path to organize correct script usage
        CWD=os.getcwd()
        
        ## create new folder for trained model and copy submit file to the folder of the architecture
        os.mkdir('TL_{}/NewModel'.format(t))        
        shutil.copy('{}/SlurmScript.sl'.format(ScriptFolder_Path),'{}/TrainModel.sl'.format(Model_Path))
        TrainJobFile=open('{}/TrainModel.sl'.format(Model_Path),'a')


        ###### Write all different tasks in the job script
        ### slicing 
        TrainJobFile.write('{}/slice_db.py -i {}/Cleaned_OverAH.smi -u {}/TL_{}/NewModel/Sliced.hr.smi -s hr -f {}/{}'.format(Model_Path,CWD,CWD,t,Model_Path,Cond_File))
        TrainJobFile.write('\n')
        
        ### training set expansion by randomizing the entries
        os.mkdir('TL_{}/NewModel/training'.format(t)) 
        TrainJobFile.write('{}/create_randomized_smiles.py -i {}/TL_{}/NewModel/Sliced.hr.smi -o {}/TL_{}/NewModel/training/ -n 50 -d multi'.format(Model_Path,CWD,t,CWD,t))
        TrainJobFile.write('\n')
        
        ### copy initial/pretrained model
        os.mkdir('TL_{}/NewModel/models'.format(t)) 
        if t == 1:
            ## for t == 1 no model of the previous iteration can be used -- because no model in t=0 was trained
            shutil.copy('{}/trained_models/{}'.format(Model_Path,PriorModel),'TL_{}/NewModel/models/model.init'.format(t))
        else:
            ## copy the last trained model from previous TL iteration
            shutil.copy('TL_{}/NewModel/models/model.trained.50'.format(t-1),'TL_{}/NewModel/models/model.init'.format(t))
    
        ### training
        TrainJobFile.write('{}/train_model.py -i {}/TL_{}/NewModel/models/model.init -o {}/TL_{}/NewModel/models/model.trained -s {}/TL_{}/NewModel/training/ -e {} -b {}'.format(Model_Path,CWD,t,CWD,t,CWD,t,epochs,batchsize))
        TrainJobFile.write('\n')
        TrainJobFile.close()
    
        ###### Submitting the Job File to the Cluster
        os.chdir('{}'.format(Model_Path))
        os.system('sbatch TrainModel.sl')
        time.sleep(5)
        TJobNames=HPC_Interaction.AllFinished($USERNAME)  ### add your hpc username here
        os.remove('TrainModel.sl')
        HPC_Interaction.MoveTJOBFiles(TJobNames,'{}/TL_{}/NewModel/'.format(CWD,t),'Training')
        os.chdir(CWD)
