# SpiroTransferLearning
Data and scripts for the publication "Machine Learning Driven Design of Spiropyran Photoswitches" - see [PlaceHolderLink]

The files consist of the full workflow including all scripts as well as all sampled and evaluated moleculs of the three individual transfer learning runs. 
The baseline model and the three fine-tuned models are available, too. Part of the workflow is the REINVENT scaffold-decorator which can be found here [github: reinvent-scaffold-decorator](https://github.com/undeadpixel/reinvent-scaffold-decorator)

________

## Data: 
Each of the files in the ```/Data/``` folder has a list of sampled molecules as a .xyz file with SMILES strings and our design criteria (addressability and thermostability). The geometries were optimized with xTB following the procedure described in the paper. For each of the molecules there is a open (spiropyran) and closed (merocyanin) form. The name of the file indicate if the structures are open or closed (```ClosedGeoms_``` or ```OpenGeoms_``` prefix), the name ```_RunX_``` indicated which of the three individual runs the molecules come from and the filename ending (```_TL0``` or ```_TL1```) indicate if the molecules are sampled from the baseline model (TL0) or the fine-tuned model (TL1).

________

## Model: 
The trained baseline and fine-tuned model(s) for the three individual runs can be found on Zenodo.org [link](https://zenodo.org/records/14011804)

________
