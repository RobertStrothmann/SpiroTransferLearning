# SpiroTransferLearning
Data and scripts for the publication "Machine Learning Driven Design of Spiropyran Photoswitches" - see [PlaceHolderLink]

The files consist of the full workflow including all scripts as well as all sampled and evaluated moleculs of the three individual transfer learning runs. 
The baseline model and the three fine-tuned models are available, too. Part of the workflow is the REINVENT scaffold-decorator which can be found here [github: reinvent-scaffold-decorator](https://github.com/undeadpixel/reinvent-scaffold-decorator)

________

## Data: 
Each of the files in the ```/Data/``` folder has a list of sampled molecules as a .xyz file with SMILES strings and our design criteria (addressability and thermostability). The geometries were optimized with xTB following the procedure described in the paper. For each of the molecules there is a open (spiropyran) and closed (merocyanin) form. The name of the file indicate if the structures are open or closed (```ClosedGeoms_``` or ```OpenGeoms_``` prefix), the name ```_RunX_``` indicated which of the three individual runs the molecules come from and the filename ending (```_TL0``` or ```_TL1```) indicate if the molecules are sampled from the baseline model (TL0) or the fine-tuned model (TL1).

________

## Model: 
The trained baseline and fine-tuned model(s) for the three individual runs can be found on [Zenodo.org](https://zenodo.org/records/14011804)

________

## Workflow: 
The full workflow is in the python script ```Run_TransferLearning.py```. The script consist of *7 major steps* that are repeated within one transfer learning iteration. Since the script is tailored for a specific use case most of the I/O parts and HPC related parts are moved out of the script. If this is interessting to you, you can read the last part of this paragraph (HPC usage). In the following a short summary of all the subroutines are listed: 

- Step 1 - Training the REINVENT Scaffold decorator: \n
  In this step the model is trained following a REINVENT scaffold- decorator training script. The ```PerformTraining()``` function in the ```Training_Part.py``` script performs a series of REINVENT scripts, namely      slicing of the training set, creating the randomized training set and training the model. 


- Step 2 - Sampling from a trained model:

- Step 3: Create .xyz structures
  
- Step 4: Optimize with xTB & lowest energy conformer assignment
  
- Step 5: Calculate excited states via TDDFT
  
- Step 6: Calculate design criteria (addressability and thermostability)
  
- Step 7: Define pareto points based on a non-dominated sorting algorithm


### HPC usage:
