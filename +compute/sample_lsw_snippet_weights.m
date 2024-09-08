function sample_lsw_snippet_weights(recompute,sampleNumber,Parray,Narray,permutations,theta)
%SAMPLE_LSW_AVERAGES this method computes Eq. 11 from the
%manuscript, i.e., the ergodic average of a_p(x), for each orbit p. In
%particular, this method computes a_p(x) for every collection of orbits
%specified by Parray, and for every integral duration specified by Narray.
%This method computes these values over a single chaotic trajectory,
%labelled by sampleNumber.
%
% Inputs:
%
%   sampleNumber    : integer describing which chaotic trajectory to
%                   compute averages over
%   theta           : value of variance to use for the Guassian kernal
%   Pmax            : the largest library size being used
%   Narray          : array of durations uses to compute the integral in
%                   time, in units of time steps. Narray must be
%                   monotonically increasing.

filename = sprintf('localdata/snippets/lsw/weights%g.mat',sampleNumber);
if isfile(filename)&&~recompute

    fprintf('snippet LSW weights for sample %g already exists. \n',sampleNumber);

else

    fprintf('computing snippet lsw averages from the sample number %g...\n',sampleNumber)

    % load in chaotic sample trajectory
    sample = load(['./localdata/chaos/sample',num2str(sampleNumber),'.mat'],'x','y','z');
    sample = [sample.x' sample.y' sample.z'];

    % allocate array to store averages in
    P = max(Parray);
    averages = zeros(P,numel(Narray));

    % for every p = 1,...,max(Parray)
    str = '';
    for p = 1:P

        % load in orbit p
        snippet = load(sprintf('localdata/snippets/snippet%g.mat',p));

        D = pdist2([snippet.x snippet.y snippet.z], sample);
        ap = trapz(exp(-1/(4*theta).*D.^2),1)/(numel(snippet.x)-1);

        ap = cumsum(ap);
        averages(p,:) = ap(Narray)./Narray;

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',p,P);
        fprintf(str);

    end

    fprintf('computing snippet lsw weights from the sample number %g...\n',sampleNumber)

    load('localdata/snippets/lsw/correlations.mat','K');
    R = size(permutations,2);
    N = numel(Narray);
    w = cell(P,1);
    ridgeparameter = 1e-6;
    str = '';
    for i = 1:P

        p = Parray(i);
        I = eye(p);
        w{i} = zeros(p,R,N);

        for r = 1:R

            ind = permutations(1:p,r);
            kernel = K(ind,ind);
            w{i}(:,r,:) = (kernel+ridgeparameter*I) \ averages(ind,:);

        end

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',i,P);
        fprintf(str);

    end

    % save out data
    save(filename,'w','theta');
    fprintf('saved results to `%s`\n',filename)

end

end


