# MM_bOED_HP-NMR
Code developed for the paper "A Bayesian Study to Optimally Unveil LDH Kinetic Mechanisms In Vivo". 

## Directories: ##


## Scripts: ##

-  ***GenerateAverageData***: Notebook to generate average data from experiment traces for the experiments considering only 3.2 mM pyruvate with changing conditions (number of cells, membrane integrity and CO2 availability). The code also extrapolates the HP-Pyruvate concentration at t=0 (mixing of pyruvate with cells) from its T1 decay value.
-  ***GenerateAverageData_DiffPyr***: Notebook to generate average data from experiment traces for the experiments considering 3 million cells with changing pyruvate concentrations. The code also extrapolates the HP-Pyruvate concentration at t=0 (mixing of pyruvate with cells) from its T1 decay value.
- ***PlotRawData***: Notebook to plot raw data traces for experiments. 
- ***VisualiseControlsHP***: Notebook to extract and plot the mean and standard deviation for the control experiments (0 mM pyruvate and 0 cells at 3.2 mM pyruvate).
- ***EntropyEstimationPosterior***: Julia script with the functions to estimate Entropy from multivariate priors/posterior distributions. See "Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit". 
