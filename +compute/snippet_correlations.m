function snippet_correlations(recompute,theta,Pmax)
%SNIPPET_CORRELATIONS This method computes the measure correlation
%matrix, K_{pq}, of snippets. 

filename = './localdata/snippets/lsw/correlations.mat';
if isfile(filename)&&~recompute

    obj = load(filename);
    if obj.theta == theta && size(obj.K,1)==Pmax
        fprintf('snippet correlations already exist, skipping... \n');
        return
    end
end

K = zeros(Pmax);
fprintf('computing correlations between %g snippets...\n', Pmax)
str = '';
for p = 1:Pmax
   
    % load in data of orbit p
    snippetp =load(sprintf('localdata/snippets/snippet%g.mat',p));
    snippetp = [snippetp.x snippetp.y snippetp.z];

    parfor q = p:Pmax

        % load in data of orbit p
        snippetq = load(sprintf('localdata/snippets/snippet%g.mat',q));
        snippetq = [snippetq.x snippetq.y snippetq.z];

        % compute the integral of the Gaussian kernel 
        distance = pdist2(snippetp, snippetq);
        G = exp(-1/(4*theta).*distance.^2);
        K(p,q) = compute.snippet_mean(compute.snippet_mean(G,2),1);

    end

    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',p,Pmax);
    fprintf(str);

end

% the above calculation only computed the lower triangular elements of
% K_pq. Since the object is symmetric, lets copy the elements over to
% the upper triangular elements as well.
K = K+K'-diag(diag(K));

% save data out
save(filename,'K','theta');
fprintf('saved results to `%s`\n',filename)

end