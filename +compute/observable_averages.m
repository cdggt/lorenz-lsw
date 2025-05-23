function observable_averages(recompute,observables,Pmax, Smax)
%OBSERVABLE_AVERAGES this methods computes all finite-time chaotic
%averages, snippet averages, and periodic orbit averages for each
%observable in the observables cell array. It also appends the lyapunov
%exponent over chaos and orbits to the results. For snippets, the lyanpunov
%exponent is set to nan. 

filename = 'localdata/predictions/averages.mat';
if isfile(filename)&&~recompute

    fprintf('averages have already been agregated, skipping... \n');

else

    fprintf('computing observable averages...\n')

    %% compute the average of each observable over each orbit and each snippet

    nObs = numel(observables);
    orbit_obs_averages = zeros(Pmax,nObs+2);
    snippet_obs_averages = zeros(Pmax,nObs+2);

    str = '';
    for p = 1:Pmax


        orbit = load(sprintf('data/orbits/orbit%g.mat',p),'x','y','z','floquetexponent');
        for o = 1:nObs
            obs = observables{o}(orbit.x,orbit.y,orbit.z);
            orbit_obs_averages(p,o) = compute.orbit_mean(obs);
        end
        lambda = orbit.floquetexponent;
        % add unstable lyapunov exponent as an extra column
        orbit_obs_averages(p,nObs+1) = lambda;
        % add stable lyapunov exponent as an extra column
        orbit_obs_averages(p,nObs+2) = -(lambda+10+8/3+1);
        % add Kaplan-Yorke dimension as an extra column
        orbit_obs_averages(p,nObs+3) = nan;

        snippet = load(sprintf('localdata/snippets/snippet%g.mat',p),'x','y','z','period');
        for o = 1:nObs
            obs = observables{o}(snippet.x,snippet.y,snippet.z);
            snippet_obs_averages(p,o) = compute.snippet_mean(obs,1);
        end
        % add nan as an extra column for unstable lyapunov exponent
        snippet_obs_averages(p,nObs+1) = nan;
        % add nan as an extra column for stable lyapunov exponent
        snippet_obs_averages(p,nObs+2) = nan;
        % add nan as an extra column for Kaplan-Yorke dimension
        snippet_obs_averages(p,nObs+3) = nan;

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',p,Pmax);
        fprintf(str);

    end

    %% compute the average of each observable over each chaotic sample

    sample_obs_averages = zeros(Smax,nObs+2);
    sample_obs_variances = zeros(Smax,nObs+2);
    str = '';
    for s = 1:Smax

        sample = load(sprintf('localdata/chaos/sample%g.mat',s),'x','y','z');
        for o = 1:nObs
            obs = observables{o}(sample.x,sample.y,sample.z);
            sample_obs_averages(s,o) = mean(obs);
            sample_obs_variances(s,o) = var(obs);
        end
        % add unstable lyapunov exponent as an extra column
        sample_obs_averages(s,nObs+1) = 0.90566;
        sample_obs_variances(s,nObs+1) = 0;
        % add stable lyapunov exponent as an extra column
        sample_obs_averages(s,nObs+2) = -14.57233;
        sample_obs_variances(s,nObs+2) = 0;
        % add stable lyapunov exponent as an extra column
        sample_obs_averages(s,nObs+3) = 2.062149;
        sample_obs_variances(s,nObs+3) = 0;

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',s,Smax);
        fprintf(str);

    end

    % save out data
    save(filename,'orbit_obs_averages','snippet_obs_averages','sample_obs_averages','sample_obs_variances');
    fprintf('saved results to `%s`\n',filename)


end

end