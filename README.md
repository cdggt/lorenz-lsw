# README of repository `lorenz-lsw`

This repository is an implments the work described in "Computing Chaotic Time-Averages from Few Periodic or Non-Periodic Orbits" by Joshua L. Pughe-Sanford, Sam Quinn, Teodor Balabanski, and Roman O. Grigoriev. It includes all scripts necessary to fully reproduce the data and figures used in the manuscript.

All periodic orbits and their stabilities were computed by Viswanath; data is available at https://websites.umich.edu/~divakar/lorenz/index.html.

## running the code

The code is run by executing `main.m`. Once this script finishes, it will plot the prediction error, Figure 2 from the paper, for each method. 

This code computes the error of each weighting scheme for all library sizes `P`, permutations `R`, chaotic trajectories `S`, and durations `N` (as described in the paper). The user has the open to change what values of `P`, `R`, `N`, and `S` are used in the header of `main.m`. 