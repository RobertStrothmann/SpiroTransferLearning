import os
import numpy as np
import argparse

### this script actually runs the REINVENT scaffold-decorator sampling procedure with a given model

def RunScaffoldSampling(Model_Path,
                        ModelToUse,
                        NumOfGPUs,
                        gpu_count,
                        ScaffoldFolder,
                        NumOfRand,
                        NumOfAtt):
   
    NumOfGPUs=int(NumOfGPUs)
    gpu_count=int(gpu_count)
    NumOfRand=int(NumOfRand)
    NumOfAtt=int(NumOfAtt)

    ## get the Scaffold Names and Scaffold SMIs
    ScaffName=np.genfromtxt('{}/ScaffLabelList.txt'.format(ScaffoldFolder),dtype=str)
    ScaffSMIs=np.genfromtxt('{}/ScaffSMIList.txt'.format(ScaffoldFolder),dtype=str)

    ## define over which area of the SMI list is iterated
    LinSpace=np.round(np.linspace(0,len(ScaffSMIs),int(NumOfGPUs)))

    if gpu_count==(NumOfGPUs-1):
        Start=int(LinSpace[gpu_count])
        Stop=None
    else:
        Start=int(LinSpace[gpu_count])
        Stop=int(LinSpace[gpu_count+1])
    
    if gpu_count==0:
        Start=int(LinSpace[0])
        Stop=int(LinSpace[1]) 
    
    ScaffSMIs_SubSpace=ScaffSMIs[Start:Stop]
    ScaffName_SubSpace=ScaffName[Start:Stop]
    
    ## loop over all possible SMIs and submit job
    for i in range(len(ScaffSMIs_SubSpace)):
        writing_file = open("scaffold_{}.smi".format(gpu_count), "w")
        writing_file.write(ScaffSMIs_SubSpace[i])
        writing_file.close()

        Name=ScaffName_SubSpace[i]
        os.system('spark-submit --driver-memory=125g {}/sample_scaffolds.py -m {} -i scaffold_{}.smi -o {}.csv --output-format=csv -r {} -n {} -d multi'.format(Model_Path,ModelToUse,gpu_count,Name,NumOfRand,NumOfAtt))



# parsing arguments 
parser = argparse.ArgumentParser(description='Samples a List of Scaffolds according to RNN Scaffold Decoration based on a Model Path and default 1024/1024 setups.')

parser.add_argument("-i", dest="Model_Path", required=True,
    help="Path to .model file", metavar="Path to model folder")

parser.add_argument("-m", dest="ModelToUse", required=True,
    help="Model that is used in sampling", metavar="Model that is used in sampling")

parser.add_argument("-gpu", dest="GPUNum", required=True,
    help="Number of GPUs available", metavar="Number of GPUs available")

parser.add_argument("-num_gpu", dest="gpu_count", required=True,
    help="GPU number for iteration", metavar="GPU number for iteration")

parser.add_argument("-scaffFolder", dest="Scaff_Folder", required=True,
    help="Path to scaffold folder", metavar="Path to scaffold folder")

parser.add_argument("-r", dest="NumOfRand", required=True,
    help="Number of randomization steps", metavar="Number of randomization steps")

parser.add_argument("-n", dest="NumOfAtt", required=True,
    help="Number of attempts", metavar="Number of attempts")

args = parser.parse_args()


## program run
RunScaffoldSampling(args.Model_Path,
                    args.ModelToUse,
                    args.GPUNum,
                    args.gpu_count,
                    args.Scaff_Folder,
                    args.NumOfRand,
                    args.NumOfAtt)
