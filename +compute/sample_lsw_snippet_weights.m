function sample_lsw_snippet_weights(recompute,sampleNumber,Parray,Narray,permutations,theta)
%SAMPLE_LSW_SNIPPET_WEIGHTS this method computes the LSW weights of snippets,
%for each P, R, and N, at a specific sample index S=sampleIndex. 

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
        ap = compute.snippet_mean(exp(-1/(4*theta).*D.^2),1);

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


