# SURE_OMT
This repository contains the code to reproduce all the figures from 
*Döhler, S., Meah, I., and Roquain, E. (2021). Online multiple testing with super-uniformity reward*.
 
## Organization of the repository
A package with the methods introduced in the paper is under construction. 
The current version of the package (used to get the numerical results in the paper) is provided here,
in OnlineSuperUnif/ folder. You must install the package to reproduce the experiments. 
```r
# install.packages("devtools") # If devtools not installed
devtools::install("OnlineSuperUnif")
```
The main/ folder contains the scripts for the experiments,
either on simulated data or on real data (see the README in it).
Parameters of the experiments are set using .json files contained in the config/ folder. 
The current parameters are the ones used for the figures provided in the paper.
Running the experiments in the main/ folder will provide Figures 6, to 13, 16, and 17 in the figures/ folder and the associated data in the data/ folder.
To launch the experiments type in a terminal
``` 
bash launch_fwerxp.sh
bash launch_mfdrxp.sh
```
Finally, in the figures/ folder and the corresponding folder in it, we provide Rmds to reproduce Figures 1, 2, 3, 4, 14, and 15. 

