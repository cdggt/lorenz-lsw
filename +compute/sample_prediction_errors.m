function sample_prediction_errors(recompute,sampleIndex,Parray,Narray)
%SAMPLE_PREDICTIONS_ERRORS this method computes takes the weights, which 
% have been saves out, and computes E_rel for each weighting scheme, over 
% each value of P, R, and N, at a specific value of S=sampleIndex. 

filename = sprintf('./localdata/predictions/errors%g.mat',sampleIndex);

if isfile(filename)&&~recompute

    fprintf('errors have already been agregated, skipping... \n');

else

    fprintf('computing errors...\n')
    load('data/library_permutations.mat','permutations');
    load('localdata/predictions/averages.mat','orbit_obs_averages','snippet_obs_averages','sample_obs_averages','sample_obs_variances');
    % for each observable, get aggegrate average and std, over all S
    sample_means = mean(sample_obs_averages,1)';
    sample_stds = sqrt(mean(sample_obs_variances,1))';

    % make sure that means that are analytically zero are actually zero
    sample_means(mean(sample_obs_averages,1)<(10^6)^(-1/2))=0;
    
    % set observables with no variance to have a variance of 1
    sample_stds(sample_stds==0)=1;

    %% allocate memory for arrays
    nObs = numel(sample_means);
    N = numel(Narray);
    R = size(permutations,2);
    P = numel(Parray);

    orbit_lsw_tikhonov_error= nan(P,R,N,nObs);
    orbit_lsw_convex1_error = nan(P,R,N,nObs);
    orbit_lsw_convex2_error = nan(P,R,N,nObs);
    orbit_markov_error      = nan(P,R,N,nObs);
    orbit_uniform_error     = nan(P,R,nObs);
    orbit_pot_error         = nan(P,R,nObs);

    snippet_lsw_tikhonov_error= nan(P,R,N,nObs);
    snippet_lsw_convex1_error = nan(P,R,N,nObs);
    snippet_lsw_convex2_error = nan(P,R,N,nObs);
    snippet_markov_error      = nan(P,R,N,nObs);
    snippet_uniform_error     = nan(P,R,nObs);

    %% load in weights

    orbit_pot_weights = load(sprintf('./localdata/orbits/pot/weights.mat'));
    orbit_pot_weights = orbit_pot_weights.w;
    orbit_markov_weights = load(sprintf('./localdata/orbits/markov/weights%g.mat',sampleIndex));
    orbit_markov_weights = orbit_markov_weights.w;
    orbit_lsw_weights = load(sprintf('./localdata/orbits/lsw/weights%g.mat',sampleIndex));

    snippet_markov_weights = load(sprintf('./localdata/snippets/markov/weights%g.mat',sampleIndex));
    snippet_markov_weights = snippet_markov_weights.w;
    snippet_lsw_weights = load(sprintf('./localdata/snippets/lsw/weights%g.mat',sampleIndex));

    %% compute errors

    % for every library size P
    str = '';
    for j = 1:numel(Parray)

        p = Parray(j);

        % for every permutation
        for r = 1:R

            % draw a random permutation
            ind = permutations(1:p,r);

            % pot prediction error 
            weights = orbit_pot_weights{j}(:,r);
            predictions = orbit_obs_averages(ind,:)'*weights;
            orbit_pot_error(p,r,:) = abs(predictions-sample_means)./sample_stds;

            % uniform prediction error 
            predictions = mean(orbit_obs_averages(ind,:),1)';
            orbit_uniform_error(p,r,:) = abs(predictions-sample_means)./sample_stds;
            predictions = mean(snippet_obs_averages(ind,:),1)';
            snippet_uniform_error(p,r,:) = abs(predictions-sample_means)./sample_stds;

            % for every sample duration
            for n = 1:N

                % markov prediction error
                weights = orbit_markov_weights{j}(:,r,n);
                predictions = orbit_obs_averages(ind,:)'*weights;
                orbit_markov_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);
                weights = snippet_markov_weights{j}(:,r,n);
                predictions = snippet_obs_averages(ind,:)'*weights;
                snippet_markov_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);

                % lsw prediction error
                weights = orbit_lsw_weights.w_tikhonov{j}(:,r,n);
                predictions = orbit_obs_averages(ind,:)'*weights;
                orbit_lsw_tikhonov_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);
                weights = orbit_lsw_weights.w_convex1{j}(:,r,n);
                predictions = orbit_obs_averages(ind,:)'*weights;
                orbit_lsw_convex1_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);
                weights = orbit_lsw_weights.w_convex2{j}(:,r,n);
                predictions = orbit_obs_averages(ind,:)'*weights;
                orbit_lsw_convex2_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);

                weights = snippet_lsw_weights.w_tikhonov{j}(:,r,n);
                predictions = snippet_obs_averages(ind,:)'*weights;
                snippet_lsw_tikhonov_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);
                weights = snippet_lsw_weights.w_convex1{j}(:,r,n);
                predictions = snippet_obs_averages(ind,:)'*weights;
                snippet_lsw_convex1_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);
                weights = snippet_lsw_weights.w_convex2{j}(:,r,n);
                predictions = snippet_obs_averages(ind,:)'*weights;
                snippet_lsw_convex2_error(p,r,n,:) = (abs(predictions-sample_means)./sample_stds);

            end

        end

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',p,P);
        fprintf(str);

    end

    %% save out results so that they may be plotted
    filename = sprintf('./localdata/predictions/errors%g.mat',sampleIndex);
    save(filename,'orbit_lsw_tikhonov_error','orbit_lsw_convex1_error','orbit_lsw_convex2_error','orbit_markov_error','orbit_uniform_error','orbit_pot_error','snippet_lsw_tikhonov_error','snippet_lsw_convex1_error','snippet_lsw_convex2_error','snippet_markov_error','snippet_uniform_error');
    fprintf('saved results to `%s`\n',filename)

end

end
