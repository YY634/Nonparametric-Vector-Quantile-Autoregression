
title: Nonparametric Vector Quantile Autoregression

**Note**: Please set the provided "reproducibility_materials" as the working directory. 
Under the working directory  (`reproducibility_materials`), there should be two sub-directories `simulations` and `real_data_analysis`. 
Please also download the folder `preprocessed_data` from the GitHub page 
https://github.com/YY634/Nonparametric-Vector-Quantile-Autoregression/tree/main/real_data_analysis/preprocessed_data
and place it under
`reproducibility_materials/real_data_analysis/`
so that the resulting path is `reproducibility_materials/real_data_analysis/preprocessed_data/`.
These are the preprocessed EEG recordings (for the real data analysis in Section 4.2). 
Due to their large size and to avoid license (data privacy) issue, we cannot include these datasets directly in the reproducibility package. Instead, we have made them available through the GitHub page provided above.  

Now the folder `reproducibility_materials` contains:
1. The R codes to conduct simulations, which replicate the simulation results.
2. The R codes for real data analysis via implementing nonparametric VQAR.
3. Preprocessed real datasets that are ready to use, and the code to preprocess the original EEG recordings.  



The sub-directory `reproducibility_materials/simulations/` contains all the R scripts needed in Section 4.1 and Appendix B.1 of the paper. 
The procedure of simulations is described in Workflow.Rmd and Workflow_simulation.txt. 

The sub-directory `reproducibility_materials/real_data_analysis/` contains all the R scripts and datasets needed in Section 4.2 of the paper. 
The procedure of real data analysis is described in Workflow.Rmd and Workflow_realdata.txt.

The purpose of each R script has been explained in the script itself and in Workflow.Rmd, Workflow_simulation.txt and Workflow_realdata.txt. 




