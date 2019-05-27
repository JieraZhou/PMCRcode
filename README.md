# PMCRcode
Precision medicine with competing risks outcome
These R codes are part of supplementary documents for the manuscript "On Restricted Optimal Treatment Regime Estimation for Competing Risks Data" by Jie Zhou, Jiajia Zhang and Wenbin Lu
This project contains 5 R code files and a pseudo HIV dataset. All files need to be downloaded and saved locally and the directory of the location need to be updated in the code respectively.

The file PMCRfun.R contains functiones needed to reproduce the simulation and fit the pseudo real data.
The simulation.R file containes the code for reproducing the simulation results in the manuscript, the read_fun.R has functions used to read the simulation results, and read.R is used to get the summary of the simulation results.
The pseudoHIV.txt is a fake dataset generated based on the HIV real data in the manuscript, and the pseudoHIVfit.R is the code used to get the estimated unrestricted and restricted treatment regimes for the pseudo HIV data.
