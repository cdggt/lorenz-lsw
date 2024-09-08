function chaotic_sample(recompute,sampleNumber,seed,N)
%CHAOTIC_SAMPLE this method computes a chaotic trajectory of
%the lorenz attractor. To make sure that everything in our manuscript is
%reproducible, we computed all 256 sample trajectories at set seeds of the
%matlab(r2023a)'s psuedorandom number generator.
%
% Inputs:
%
%   sampleNumber    : integer describing which chaotic trajectory to
%                   compute
%   seed            : seed for the random number generator
%   N               : number of samples to save out per trajectory

% only compute trajectory if it does not already exist
filename = ['./localdata/chaos/',sprintf('sample%g',sampleNumber),'.mat'];
if isfile(filename)&&~recompute

    fprintf('sample number %g already exists. \n',sampleNumber);

else
    % set the seed. This is to ensure the data created here matches the data in the paper exactly.
    rng(seed+sampleNumber);

    % set data size parameters

    % time step of the RK4 integrator
    timestep = 2e-3;

    % before collecting data, integrate for this many time units to make sure
    % the chaotic states lies on the chaotic attractor
    settlingTime = 25;

    % save out the state of the chaotic trajectory approximately every 'dtPerDelta*dt' time units
    dtPerReadout = 1000;

    %% Collect Data

    fprintf('computing chaotic sample number %g\n',sampleNumber)

    % allocate arrays
    x = zeros(3,N);
    t = zeros(N,1);

    % initialize on a random point. Integrate for the bit tho to make sure
    % that this point settles onto the chaotic attractor before taking data
    fprintf('\t integrating onto the attractor...\n');
    x(:,1) = [5*rand; 5*rand; 20];
    for i = 1:ceil(settlingTime/timestep)
        x(:,1) = rk4(x(:,1),timestep);
    end

    % actually take data
    tic % begin timer
    runtime = 0;
    init = false;
    for n = 2:N

        x(:,n) = x(:,n-1);
        for k = 1:dtPerReadout
            x(:,n) = rk4(x(:,n),timestep);
        end
        t(n) = t(n-1)+dtPerReadout*timestep;

        % print at most once every 60 seconds
        if toc > 60 || init
            runtime = runtime+toc;
            percentdone = n/N*100;
            timeleft = (100-percentdone)/percentdone*runtime;
            fprintf('\t collecting data, %g%% complete (%g min left) \n',round(percentdone), max(timeleft,0));
            init=false;
            tic;
        end

    end

    % save out chaotic trajectory
    z = x(3,:);
    y = x(2,:);
    x = x(1,:);
    save(filename,'x','y','z','t');
    fprintf('saved results to `%s`\n',filename)

end

end
