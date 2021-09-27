# SURE_OMT
This repositery contains the code to reproduce all the figures from 
*Döhler, S., Meah, I., and Roquain, E. (2021). Online multiple testing with super-uniformity reward*.
 
## Organization of the repository
A package, with the methods introduced in the paper, is under construction. 
The current version of the package (used to get the numerical results in the paper) is provided here,
in the folder OnlineSuperUnif/. You must, first, install the package to reproduce the experiments. 
```{r, message=FALSE, warning=FALSE, results=FALSE}
# install.packages("devtools") # If devtools not installed
devtools::install("OnlineSuperUnif")
```
The main/ folder contains the scripts for the experiments (either on simulated data or on real data).
The parameter for the experiment are set using .json file contained in the config/ folder. 
The current parameters are the one used for the figures provided in the paper.
Running the experiments in the main/ folder will provide Figures .... in the figures/ folder and the associated data in the data/ folder.
To launch the experiments type in a terminal
``` 
bash launch_fwerxp.sh
bash launch_mfdrxp.sh
```
Finally, for Figures ... we provide Rmd to reproduce them still in the figures/ folder and the corresponding folder in it. 

