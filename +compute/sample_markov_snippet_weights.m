function sample_markov_snippet_weights(recompute,sampleNumber,Parray,Narray,permutations)
%SAMPLE_MARKOV_SNIPPET_WEIGHTS this method computes the Markov weights of 
% snippets, for each P, R, and N, at a specific sample index S=sampleIndex. 

filename = sprintf('localdata/snippets/markov/weights%g.mat',sampleNumber);
if isfile(filename)&&~recompute

    fprintf('orbit markov weights for sample %g already exists. \n',sampleNumber);

else

% load in chaotic sample trajectory
sample = load(['./localdata/chaos/sample',num2str(sampleNumber),'.mat'],'x','y','z');
sample = [sample.x' sample.y' sample.z'];

fprintf('computing snippet Markov averages for sample number %g...\n',sampleNumber)

% allocate memory for arrays
Pmax = numel(Parray);
D = nan(Pmax,size(sample,1));

str = '';
for p = 1:Pmax

    snippet = load(sprintf('localdata/snippets/snippet%g.mat',p));
    D(p,:) = min(pdist2([snippet.x snippet.y snippet.z], sample),[],1);

    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',p,Pmax);
    fprintf(str);

end

fprintf('computing snippet Markov weights for sample number %g...\n',sampleNumber)

P = numel(Parray);
R = size(permutations,2);
N = numel(Narray);
w = cell(P,1);
str = '';
for i = 1:P
    p = Parray(i);
    w{i} = zeros(p,R,N);
    for r = 1:R

        ind = permutations(1:p,r);
        [~,closest] = min(D(ind,:),[],1);

        for n = 1:N
                
            w{i}(:,r,n) = hist(closest(1:Narray(n)),1:p)'/Narray(n);

        end

    end

    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',i,P);
    fprintf(str);

end

% save out data
save(filename,'w');
fprintf('saved results to `%s`\n',filename)

end

end