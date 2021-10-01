This folder contains the code to obtain the results on simulated data or real data. The scripts are separated in fwer/ and mfdr/ folder, to get the results on fwer controlling procedures and on mfdr controlling procedures (the figures/ and data/ folders are organized in the same manner).

Each script corresponds to one parameter of study: the position of the signal, its proportion (piA), its strength (p3), lambda, the number of individuals N, and the kernel bandwidth h (see the paper for more details). 

The scrips named "impc*"  are for the real data application. We use data provided by Karp et al. [here](https://zenodo.org/record/260398#.YVa84kZBzUI), but for space reason, the initial .csv file (ReproducibleCode/Figure 4/SDofGenotypeEffect_processedData_categorical.csv) could not be added to this repository. Instead, we added a sub-data frame containing only the variables that we need. At the beginning of the "impc*" scripts you can see the code used to extract these variables and create the sub-data frame. 

In each experiment run, the results are stored in csv files saved in the data/ folder and plots are created from this data and saved in the figures/ folder.