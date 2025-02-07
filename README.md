# README of repository `lorenz-lsw`

This repository implements the work described in "Computing Chaotic Time-Averages from Few Periodic or Non-Periodic Orbits" by Joshua L. Pughe-Sanford, Sam Quinn, Teodor Balabanski, and Roman O. Grigoriev. It includes all scripts necessary to fully reproduce the data and figures used in the manuscript.

All periodic orbits and their stabilities were computed by Viswanath; data is available at https://websites.umich.edu/~divakar/lorenz/index.html.

If you ever have questions about this code, or want help adapting it to your own use case, please feel free to reach out to me at jpughesanford@{gmail.com,gatech.edu}.

## Running the Code

The code is run in full simply by executing `main.m`. Once this script finishes, it will generate all figures in the paper and save the relevant data to disk.

Note that the purpose of this code base goes beyond computing LSW weights. This code quantifies the accuracy of LSW, Markov, POT, and uniform weighting schemes over four hyper-parameters, $P$, $R$, $S$, and $N$. See manuscript for more details.

These variables are defined in `main.m` via the variables:
- `Parray`: an array containing the values of $P$ to be tested. We use `Parray = 1...125`
- `R`: the number of library permutations to compute the accuracy over. We use $R=256$
- `S`: the number of chaotic trajectories to measure the accuracy over. We use $S=256$
- `Narray`: an array containing the values of $N$ to be tested. We use $Narray=[10, 10^2, ..., 10^6]$

Given the number of variables accuracy is being computed over, this main script can take a long time to run. We computed everything in parallel on a cluster. The code base, as written here, runs in series on a single process. If you would like to run this code on your local machine, without having to wait too long, we suggest making the values of `S`, `R`, or `P` smaller. Alternatively, note that all for loops over `S` can be run in parallel across multiple processes using Matlab's `parfor` loop. 

## Folder Structure

Once `main` has finished running, the project folder will be populated with the following folders 

- `/media` - figures will be saved out here
- `/data/orbits/` - Viswanath's orbits `orbit{index}.mat` are saved here, where $\text{index} = 1,\dots,125$. 
- `/localdata/chaos/` - chaotic trajectories `sample{index}.mat` will be populated here, where $\text{index} = 1,\dots,S$. 
- `/localdata/orbits/lsw` - data related to the LSW weights of orbits will be populated here
- `/localdata/orbits/markov` - data related to the Markov weights of orbits will be populated here
- `/localdata/orbits/pot` - data related to the LSW weights of orbits will be populated here
- `/localdata/snippets/` - all of the snippets in the snippet library will be populated here
- `/localdata/snippets/lsw` - data related to the LSW weights of snippets will be populated here
- `/localdata/snippets/markov` - data related to the Markov weights of snippets will be populated here
- `/localdata/predictions` - data related to computing method errors will be populated here

## Data Structures

Let lowercase $p$, $r$, $n$, and $s$ denote specific values $p\in$ `Parray`, $n\in$ `Narray`, $r\in [1,2, \dots, R]$, and $s\in [1,2, \dots, S]$. Also, define 
- $N_0$: the number of elements in `Narray`
- $N_p$: the number of snapshots in time saved out from orbit $p$
- $N_\text{max}$: the largest value of $N$ in `Narray`. 
- $P_0$: the number of elements in `Parray`
- $P_\text{max}$: the total number of known orbits or snippets. Here, $P\text{max}=125$. 
- $B_0$: the number of observables in $\mathcal{B}$. We use $\mathcal{B} = \{1,x,y,z,x^2,xy,xz,y^2,yz,z^2\}$ such that $B_0=10$.

### trajectory data 

Every chaotic trajectory .mat file in `/localdata/chaos/` has the following fields:
- `x`: an $[N_\text{max} \times 1]$ array of the $x$-coordinate, over the chaotic trajectory
- `y`: an $[N_\text{max} \times 1]$ array of the $y$-coordinate, over the chaotic trajectory
- `z`: an $[N_\text{max} \times 1]$ array of the $z$-coordinate, over the chaotic trajectory
- `t`: an $[1 \times N_\text{max}]$ array of the $t$-coordinate, over the chaotic trajectory
    
Every periodic orbit .mat file in `/data/orbits/` has the following fields:
- `x`: an $[N_p \times 1]$ array of the $x$-coordinate, over the periodic orbit
- `y`: an $[N_p \times 1]$ array of the $y$-coordinate, over the periodic orbit
- `z`: an $[N_p \times 1]$ array of the $z$-coordinate, over the periodic orbit
- `period`: the period of the orbit, $T_p$
- `floquetexponent`: the top Floquet exponent of the orbit, $\lambda_p$
- `topologicalperiod`: the symbol length of the orbit

Every snippet .mat file in `/localdata/snippets/` has the following fields:
- `x`: an $[N_p \times 1]$ array of the $x$-coordinate, over the snippet
- `y`: an $[N_p \times 1]$ array of the $y$-coordinate, over the snippet
- `z`: an $[N_p \times 1]$ array of the $z$-coordinate, over the snippet
- `period`: the period of the snippet, $T_p$

### POT Data 

The `/localdata/orbits/pot/` folder contains a file `weights.mat`. This file  has fields:
- `w`: a $\{P_0\times 1\}$ cell array containing periodic orbit weights computed at each $p$

Each element of the `w` cell array is a $[p\times R]$ matrix. For a given $p$ and $r$, the POT weights are
>pot_weights = w{p}(:,r)

Since POT weights do not depend on the chaotic trajectory, the weights do not vary with $n$ or $s$. 

### Markov Data 

The `/localdata/{orbits,snippets}/markov/` folders each have numerous `weights{index}.mat` files, one for each value of $s$. Each file pertains to a specific chaotic trajectory. Each file contains the following fields :
- `w`: a $\{\#_P\times 1\}$ cell array containing the Markov weights for each library size $P$ in `Parray`

Each element of the `w` cell array is a $[p\times R\times N_0]$ matrix. The Markov weights computed for a given $p$, $r$, $n$, and $s$ are 
>markov_weights = w{p}(:,r,n)

where `w` is loaded from `weights{s}.mat`.

### LSW Data

The `/localdata/{orbits,snippets}/lsw/` folders each have a `correlations.mat` file containing:
- `K`: a $[P_\text{max}\times P_\text{max}]$ matrix of orbit (or snippet) correlations. This matrix is called $A_{pq}$ in the paper. 
- `theta`: the Gaussian kernel variance used to compute the correlations

The `/localdata/{orbits,snippets}/lsw/` folders each have numerous `weights{index}.mat` files, one for each value of $s$. Each file pertains to a specific chaotic trajectory. Each file contains the following fields :
- `theta`: the Gaussian kernel variance used to compute the LSW weights
- `w_tikhonov`: a $\{P_\text{max}\times 1\}$ cell array containing the LSW weights, computed using Tikhonov regularization, for each library size $P$
- `w_convex1`: a $\{P_\text{max}\times 1\}$ cell array containing the LSW weights, computed using Matlab's $\verb|lsqnonneg|$, for each library size $P$
- `w_convex2`: a $\{P_\text{max}\times 1\}$ cell array containing the LSW weights, computed using Matlab's $\verb|fmincon|$, for each library size $P$

Each element of the weight cell arrays is a $[p\times R\times N_0]$ matrix. The LSW weights computed for a given $p$, $r$, $n$, and $s$ are 
>lsw_weights_unconstrained = w_tikhonov{p}(:,r,n)\
>lsw_weights_convex_lsqnonneg = w_convex1{p}(:,r,n)\
>lsw_weights_convex_fmincon = w_convex2{p}(:,r,n)

where the weight cell arrays are loaded from `weights{s}.mat`.

### Accuracy

The `/localdata/predictions/` folder contains an `averages.mat` file containing:
- `orbit_obs_averages`: a $[P_\text{max}\times (B_0+1)]$ matrix of test observable averages over orbits. The first $B_0$ columns correspond to the averages over the elements of $\mathcal{B}$. The very last column is populated with every orbits largest Floquet exponent, $\lambda_p$.  
- `snippet_obs_averages`: a $[P_\text{max}\times (B_0+1)]$ matrix of test observable averages over snippets. The first $B_0$ columns correspond to the averages over the elements of $\mathcal{B}$. The very last column is set to NaN.  
- `sample_obs_averages`: a $[1\times (B_0+1)]$ array of matrix of test observable averages over *all* chaotic data (all $S$ trajectories). The very last column is populated with the top Lyapunov exponent of the Lorenz attractor. 
- `sample_obs_averages`: a $[1\times (B_0+1)]$ array of matrix of test observable variances over *all* chaotic data (all $S$ trajectories). The very last column is set to 1. 

The `/localdata/predictions/` folder contains numerous `errors{index}.mat` files, one for each value of $s$. Each file pertains to a specific chaotic trajectory. Each file contains the following fields:

- `orbit_pot_error`: the relative error of a POT weighting, over orbits
- `orbit_uniform_error`: the relative error of a uniform weighting, over orbits
- `snippet_uniform_error`: the relative error of a uniform weighting, over snippets
- `orbit_markov_error`: the relative error of a Markov weighting, over orbits
- `snippet_markov_error`: the relative error of a Markov weighting, over snippets
- `orbit_lsw_tikhonov_error`: the relative error of an LSW weighting, computed using Tikhonov regularization, over orbits
- `orbit_lsw_convex1_error`: the relative error of an LSW weighting, computed using $\verb|lsqnonneg|$, over orbits
- `orbit_lsw_convex2_error`: the relative error of an LSW weighting, computed using $\verb|fmincon|$, over orbits
- `snippet_lsw_tikhonov_error`: the relative error of an LSW weighting, computed using Tikhonov regularization, over snippets
- `snippet_lsw_convex1_error`: the relative error of an LSW weighting, computed using $\verb|lsqnonneg|$, over snippets
- `snippet_lsw_convex2_error`: the relative error of an LSW weighting, computed using $\verb|fmincon|$, over snippets

The fields `orbit_pot_error`, `orbit_uniform_error`, and `snippet_uniform_error` are all $[P_0\times R\times (B_0+1)]$ arrays. For a given $p$, $r$, $n$, and $s$, the value of $E_\text{max}$ the POT weights are
>POT_error = max(orbit_pot_error(p,r,1:end-1))
where `orbit_pot_error` is loaded from `errors{s}.mat`.

Notice that we compute the max only over the elements of $\mathcal{B}$, and not over the Lyanpunov exponent. 

The remaining fields are all $[P_0\times R\times \times N \times (B_0+1)]$ arrays. For a given $p$, $r$, $n$, and $s$, the relative error of, for example, the unconstrained LSW weights are
>LSW_error = max(orbit_lsw_tikhonov_error(p,r,n,1:end-1))
where `orbit_lsw_tikhonov_error` is loaded from `errors{s}.mat`.

