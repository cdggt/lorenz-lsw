close all
clear all

mkdir localdata
mkdir media

mkdir localdata/chaos

mkdir localdata/orbits
mkdir localdata/orbits/lsw
mkdir localdata/orbits/markov
mkdir localdata/orbits/pot

mkdir localdata/snippets
mkdir localdata/snippets/lsw
mkdir localdata/snippets/markov

mkdir localdata/predictions

recompute = false; % when false, this code will not recompute data whose files already exist in ./localdata

% NOTE: While this script computes the data on a single process, in series,
% the data for the paper was computed in parallel on a cluster. There is a
% lot being computed, so this script may take a long time to run. 

%% define parameters

% Define P,R,S,and N. see README. 

Parray = 1:125; % the library sizes to consider
Pmax = max(Parray);

S = 256; % the number of chaotic trajectories to compute
R = 256; % the number of library permutations to compute

Narray = 10.^(1:6); % the choatic trajectory durations to consider when computing weights

seed = 123; % this is the seed we used to generate the values in the paper
rng(seed);

%% generate a snippet library to match Viswanaths orbit library

compute.snippet_library(Pmax,seed);

%% pick what value to use for the Gaussian kernel variance

theta = 10^2;

%% compute chaotic trajectories, and LSW, MARKOV, POT weights 

% compute K_{pq} for LSW weighting 
compute.orbit_correlations(recompute,theta,Pmax); % orbits 
compute.snippet_correlations(recompute,theta,Pmax); % snippets

% compute periodic orbit weights
compute.pot_orbit_weights(recompute,Parray);

for sampleIndex = 1:S

    % compute a sample chaotic trajectory
    compute.chaotic_sample(recompute,sampleIndex,seed,max(Narray));

    % compute lsw weights, for this sample, over all p in Parray, n in Narray, and r = 1,...,R
    compute.sample_lsw_orbit_weights(recompute,sampleIndex,Parray,Narray,theta); % orbits
    compute.sample_lsw_snippet_weights(recompute,sampleIndex,Parray,Narray,theta); % snippets 

    % compute markove weights, for this sample, over all p in Parray, n in Narray, and r = 1,...,R
    compute.sample_markov_orbit_weights(recompute,sampleIndex,Parray,Narray); % orbits
    compute.sample_markov_snippet_weights(recompute,sampleIndex,Parray,Narray); % snippets

end

%% compute test observable averages, as well as the Lyapunov exponent

observables = {
    @(x,y,z) ones(size(x)),...
    @(x,y,z) x, ...
    @(x,y,z) y, ...
    @(x,y,z) z, ...
    @(x,y,z) x.*x, ...
    @(x,y,z) x.*y, ...
    @(x,y,z) x.*z, ...
    @(x,y,z) y.*y, ...
    @(x,y,z) y.*z, ...
    @(x,y,z) z.*z ...
    % compute.observable_averages will also append the lyapunov exponent
    % and Kaplan-Yorke dimension "observables" to this list.
};
compute.observable_averages(recompute,observables,Pmax,S)

%% compute E_max and E_rel for each test observable, over each trajectory

% compute E_rel for each observable, over each chaotic sample individually
for sampleIndex = 1:S
    compute.sample_prediction_errors(recompute,sampleIndex,Parray,Narray);
end

%% plot Figures

plotFigure1(recompute);
plotFigure2(recompute);
plotFigure3(Parray, R, S, Narray); 
plotFigure4(Parray, R, S, Narray);
plotTable1(Parray, R, S, Narray);