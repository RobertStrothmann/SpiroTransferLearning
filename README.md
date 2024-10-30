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

- Step 1 - Training the REINVENT Scaffold decorator: <br/>
  In this step the model is trained following a REINVENT scaffold- decorator training script. The ```PerformTraining()``` function in the ```Training_Part.py``` script performs a series of REINVENT scripts, namely      slicing of the training set, creating the randomized training set and training the model. In first iteration of the transfer learning no retraining is needed.


- Step 2 - Sampling from a trained model: <br/>
  This step samples molecules from a trained model by the use of a REINVENT scaffold-decorator script. The ```PerformSampling()``` function in the ```Sampling_Part.py``` script creates input files for the actual   
  sampling script ```RunSampling.py```. On has to give a so callled "Scaffold Folder" path in which the scaffolds as well as their labels are given. You can find an example for this folder above. 

- Step 3 - Create .xyz structures: <br/>
  Beside some I/O related commands this part provides the script to convert the sampled molecular SMILES from step 2 into .xyz files. To do so ```CSV2XYZ.py``` and ```SMI2XYZ.py``` are used to create a fixed number
  of most likely molecules as .xyz coordinates of the closed (SP) and open (MC) form. To account for different conformers the script creates not only one .xyz structure per molecule, but 10 by the use of the
  ```AllChem.EmbedMultipleConfs``` module in RdKit. 
  
- Step 4 - Optimize with xTB & lowest energy conformer assignment: <br/>
  First, a script (```GenxTB_Inputs.py```) creates input files for the following xTB geometry optimization which is then run and analyzed for non converged structures by the script ```MoveFailedRuns.py```. After
  that, a script (```GetMinConformer.py```) is used to select the lowest xTB energy conformer for a given structure and compress the other conformers. 

- Step 5 - Calculate excited states via TDDFT: <br/>
  The lowest xTB energy conformer from step 4 is used to run TDDFT calculations and obtain the excited states. To do so, the function ```GenSlurmFiles_TDDFT()``` from the script ```GenSlurmScripts.py``` is used to
  create HPC submission files
  
- Step 6 - Calculate design criteria (addressability and thermostability): <br/>
  The results (optimized geometries and excited states) from step 4 and step 5 are used to calculate the design criteria. Beside some I/O and HPC related commands the scripts ```OpenClosedAnalysis.py``` and
  ```ThermoStabilityAnalysis.py``` are used to evaluate candidate structures and write their addressability and thermostability values out. 
  
- Step 7 - Define pareto points based on a non-dominated sorting algorithm: <br/>
  In the last step, all calculated criteria from step 6 are used in a non-dominated sorting algorithm (applied in the script ```ParetoAnalysis.py```) to obtain a set of (close to) pareto front molecules. The hits are
  written down into a file and are used in the following transfer learning iteration for retraining.


### HPC usage:
