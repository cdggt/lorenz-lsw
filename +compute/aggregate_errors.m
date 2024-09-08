function aggregate_errors(Parray, S, Narray, permutations)

filename = 'localdata/predictions/errors.mat';
if isfile(filename)

    fprintf('errors have already been agregated, skipping... \n');

else

    fprintf('computing errors...\n')
    load('localdata/predictions/averages.mat','orbit_obs_averages','snippet_obs_averages','sample_obs_averages','sample_obs_variances');
    sample_means = mean(sample_obs_averages,1)';
    sample_means(mean(sample_obs_averages,1)<(10^6)^(-1/2))=0;
    sample_stds = sqrt(mean(sample_obs_variances,1))';
    sample_stds(sample_stds==0)=1;

    %% allocate memory for arrays
    N = numel(Narray);
    R = size(permutations,2);
    P = numel(Parray);

    orbit_lsw_error    = nan(P,R,N,S);
    orbit_markov_error = nan(P,R,N,S);
    orbit_uniform_error= nan(P,R);
    orbit_pot_error    = nan(P,R);

    snippet_lsw_error    = nan(P,R,N,S);
    snippet_markov_error = nan(P,R,N,S);
    snippet_uniform_error= nan(P,R);

    %% load in correlation matrices, K_{pq}

    Korbit = load('localdata/orbits/lsw/correlations.mat');
    Korbit = Korbit.K;
    Ksnippet = load('localdata/snippets/lsw/correlations.mat');
    Ksnippet = Ksnippet.K;
    ridgeparameter = 1e-6;

    %% compute errors

    orbit_pot_weights = load(sprintf('./localdata/orbits/pot/weights.mat'));
    orbit_pot_weights = orbit_pot_weights.w;

    % for every sample
    for s = 1:S

        % load in the precomputed metrics for sample s
        orbit_markov_weights = load(sprintf('./localdata/orbits/markov/weights%g.mat',s));
        orbit_markov_weights = orbit_markov_weights.w;
        snippet_markov_weights = load(sprintf('./localdata/snippets/markov/weights%g.mat',s));
        snippet_markov_weights = snippet_markov_weights.w;

        orbit_lsw_weights = load(sprintf('./localdata/orbits/lsw/weights%g.mat',s));
        orbit_lsw_weights = orbit_lsw_weights.w;
        snippet_lsw_weights = load(sprintf('./localdata/snippets/lsw/weights%g.mat',s));
        snippet_lsw_weights = snippet_lsw_weights.w;

        % orbit_lsw_averages = load(sprintf('./localdata/orbits/lsw/averages%g.mat',s));
        % orbit_lsw_averages = orbit_lsw_averages.averages;
        % snippet_lsw_averages = load(sprintf('./localdata/snippets/lsw/averages%g.mat',s));
        % snippet_lsw_averages = snippet_lsw_averages.averages;

        % for every library size p in Parray
        str = '';
        for j = 1:P

            p = Parray(j);

            % for every permutation
            for i = 1:R

                % draw a random permutation
                ind = permutations(1:p,i);
                
                if s==1

                    % compute pot prediction using p random orbits
                    weights = orbit_markov_weights{j}(:,i);
                    predictions = orbit_obs_averages(ind,:)'*weights;
                    orbit_pot_error(p,i) = max(abs(predictions-sample_means)./sample_stds);

                    % compute uniform predictions using p random orbits
                    predictions = mean(orbit_obs_averages(ind,:),1)';
                    orbit_uniform_error(p,i) = max(abs(predictions-sample_means)./sample_stds);
                    predictions = mean(snippet_obs_averages(ind,:),1)';
                    snippet_uniform_error(p,i) = max(abs(predictions-sample_means)./sample_stds);

                end

                % for every sample duration
                for n = 1:N

                    % compute markov predictions using p random orbits over Narray(n) snapshots of chaos
                    weights = orbit_markov_weights{j}(:,i,n);
                    predictions = orbit_obs_averages(ind,:)'*weights;
                    orbit_markov_error(p,i,n,s) = max(abs(predictions-sample_means)./sample_stds);
                    weights = snippet_markov_weights{j}(:,i,n);
                    predictions = snippet_obs_averages(ind,:)'*weights;
                    snippet_markov_error(p,i,n,s) = max(abs(predictions-sample_means)./sample_stds);

                    % compute lsw predictions using p random orbits over Narray(n) snapshots of chaos
                    weights = orbit_lsw_weights{j}(:,i,n);
                    predictions = orbit_obs_averages(ind,:)'*weights;
                    orbit_lsw_error(p,i,n,s) = max(abs(predictions-sample_means)./sample_stds);
                    weights = snippet_lsw_weights{j}(:,i,n);
                    predictions = snippet_obs_averages(ind,:)'*weights;
                    snippet_lsw_error(p,i,n,s) = max(abs(predictions-sample_means)./sample_stds);

                end

            end

        end


        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g\n',s,S);
        fprintf(str);

    end

    %% save out results so that they may be plotted
    filename = sprintf('./localdata/predictions/errors.mat');
    save(filename,'orbit_lsw_error','orbit_markov_error','orbit_uniform_error','orbit_pot_error','snippet_lsw_error','snippet_markov_error','snippet_uniform_error');
    fprintf('saved results to `%s`\n',filename)

end

end
