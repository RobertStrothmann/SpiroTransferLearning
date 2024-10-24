import os
import time
import shutil

#### dependent on the hpc system this script has to be adjusted !!!


## helper function which moves the finished job log files to specific destination and rename them
def MoveTJOBFiles(TJobNames,Destination,Name):
        
        TJobOut=[]
        TJobErr=[]
        for name in TJobNames:
            for i in os.listdir():
                if 'tjob.out.{}'.format(name) in i:
                    TJobOut.append(i)
                if 'tjob.err.{}'.format(name) in i:
                    TJobErr.append(i)
                    
        for i in range(len(TJobOut)):
            shutil.move(TJobOut[i],'{}/{}_{}.jobout'.format(Destination,Name,TJobNames[i]))
            shutil.move(TJobErr[i],'{}/{}_{}.joberr'.format(Destination,Name,TJobNames[i]))    
        
    



## helper function to check if the 
def AllFinished(ClusterUsername):
    ## get an inital file of running jobs
    os.system('squeue -u {} > RunJobList'.format(ClusterUsername))
    RunJob=open('RunJobList','r')
    RunningJobs=RunJob.readlines()
    RunJob.close()
    
    ## get ID of all running jobs
    TJobNames=[]
    for i in range(len(RunningJobs)):
        if i == 0:
            pass
        else:
            TJobNames.append(RunningJobs[i].split()[0])
        
    ## wait until all jobs are done
    while (len(RunningJobs) > 1):
        time.sleep(30)
        os.system('squeue -u {} > RunJobList'.format(ClusterUsername))
        RunJob=open('RunJobList','r')
        RunningJobs=RunJob.readlines()
        RunJob.close()
    os.remove('RunJobList')
    time.sleep(5)
    return(TJobNames)
