function observable_averages(recompute,observables,Pmax, Smax)

filename = 'localdata/predictions/averages.mat';
if isfile(filename)&&~recompute

    fprintf('averages have already been agregated, skipping... \n');

else

    fprintf('computing observable averages...\n')

    %% compute the averages of this observable over each orbit and each snippet

    nObs = numel(observables);
    orbit_obs_averages = zeros(Pmax,nObs+1);
    snippet_obs_averages = zeros(Pmax,nObs+1);

    str = '';
    for p = 1:Pmax


        orbit = load(sprintf('data/orbits/orbit%g.mat',p),'x','y','z','floquetexponent');
        for o = 1:nObs
            obs = observables{o}(orbit.x,orbit.y,orbit.z);
            orbit_obs_averages(p,o) = fftmean(obs);
        end
        orbit_obs_averages(p,nObs+1) = orbit.floquetexponent;

        snippet = load(sprintf('localdata/snippets/snippet%g.mat',p),'x','y','z','period');
        for o = 1:nObs
            obs = observables{o}(snippet.x,snippet.y,snippet.z);
            snippet_obs_averages(p,o) = trapz(obs)/(numel(obs)-1);
        end
        snippet_obs_averages(p,nObs+1) = nan;

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',p,Pmax);
        fprintf(str);

    end

    %% compute the averages of this observable over each chaotic sample

    sample_obs_averages = zeros(Smax,nObs+1);
    sample_obs_variances = zeros(Smax,nObs+1);
    str = '';
    for s = 1:Smax

        sample = load(sprintf('localdata/chaos/sample%g.mat',s),'x','y','z');
        for o = 1:nObs
            obs = observables{o}(sample.x,sample.y,sample.z);
            sample_obs_averages(s,o) = mean(obs);
            sample_obs_variances(s,o) = var(obs);
        end
        sample_obs_averages(s,nObs+1) = 0.90566;
        sample_obs_variances(s,nObs+1) = 1;

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',s,Smax);
        fprintf(str);

    end

    % save out data

    save(filename,'orbit_obs_averages','snippet_obs_averages','sample_obs_averages','sample_obs_variances');
    fprintf('saved results to `%s`\n',filename)


end

end