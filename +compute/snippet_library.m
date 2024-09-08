function snippet_library(Pmax,seed)
%SNIPPET_LIBRARY This method computes a library of short chaotic
%trajectories, i.e. snippets. For each periodic orbit reference in Parray,
%this method will save a snippet of equal length. This method is NOT
%designed to be intelligent, e.g., compute snippets that well cover the
%chaotic attractor. It computes snippets at random, to show that the Least 
% Squares Weighting can adjust to any snippet library.

% Inputs:
%
%   Pmax            : the largest library size being used
%   seed            : seed for the random number generator

% set the seed. This is to ensure the data created here matches the data in the paper exactly.
rng(seed+1e8); % the offset to the seed makes sure the chaotic traj generated here will not be the same as any sample chaotic trajectories

% set data size parameters

% time step of the RK4 integrator
integration_timestep = 2e-3;

% before collecting data, integrate for this many time units to make sure
% the chaotic states lies on the chaotic attractor
settlingTime = 25;

%% compute the total duration of all orbits

totalDuration = 0;
totalTimesteps = 0;
for p = 1:Pmax
    orbit = load(sprintf('data/orbits/orbit%g.mat',p));
    totalDuration = totalDuration+orbit.period;
    totalTimesteps = totalTimesteps+numel(orbit.x);
end

% set the timestep of the saved out snippets such that both sets have
% roughly equal number of time steps
saveout_timestep = totalDuration/totalTimesteps;

%% compute a chaotic trajectory of that same duration

N = ceil(totalDuration/integration_timestep);
dt = totalDuration/N;

% initialize on a random point. Integrate for the bit tho to make sure
% that this point settles onto the chaotic attractor before taking data
X0 = [5*rand; 5*rand; 20];
for i = 1:ceil(settlingTime/dt)
    X0 = rk4(X0,dt);
end

% integrate a snippet that is equal in length to the orbit
X = zeros(3,N);
X(:,1) = X0;
t = zeros(N,1);
for n = 1:N

    X(:,n+1) = rk4(X(:,n),dt);
    t(n+1) = t(n)+dt;

end

%% compute snippets as segments of this chaotic trajectory

skip = ceil(saveout_timestep/integration_timestep);
T = totalDuration/Pmax;

for p = 1:Pmax

    filename = ['./localdata/snippets/',sprintf('snippet%g',p),'.mat'];

    ind = 1:numel(t);
    ind = ind( t>=(p-1)*T & t<p*T );
    ind = ind( 1:skip:end );
    
    z = X(3,ind)';
    y = X(2,ind)';
    x = X(1,ind)';
    period = t(ind(end))-t(ind(1));

    save(filename,'x','y','z','period');

end

end
