% this is the seed we used to generate the values in the paper
seed = 123;
rng(seed);

% compute the library permutations, {P_r}, explicitly 
R = 256;
P = 125;
permutations = (1:P);
while size(permutations,1)<R
    permutations(end+1,:) = randperm(P);
    permutations = unique(permutations,'rows');
end
permutations = permutations';

% save it out. This will be saved with the Github repo just to make sure
% that if the rng function changes in MATLAB, the permutations used in this
% paperpersist.
save('localdata/library_permutations.mat','permutations');