function orbit_correlations(recompute,theta,Pmax)
%COMPUTE_ORBIT_CORRELATIONS This method computes the measure correlation
%matrix, K_{pq}, of orbits. 

filename = './localdata/orbits/lsw/correlations.mat';
if isfile(filename)&&~recompute

    obj = load(filename);
    if obj.theta == theta && size(obj.K,1)==Pmax
        fprintf('orbit correlations already exist, skipping... \n');
        return
    end

end

K = zeros(Pmax);
fprintf('computing correlations between %g orbits...\n', Pmax)
str = '';
for p = 1:Pmax
    
    % load in data of orbit p
    orbitp =load(sprintf('data/orbits/orbit%g.mat',p));
    orbitp = [orbitp.x orbitp.y orbitp.z];

    parfor q = p:Pmax

        % load in data of orbit p
        orbitq = load(sprintf('data/orbits/orbit%g.mat',q));
        orbitq = [orbitq.x orbitq.y orbitq.z];

        % compute the integral of the Gaussian kernel 
        distance = pdist2(orbitp, orbitq);
        G = exp(-1/(4*theta).*distance.^2);
        K(p,q) = compute.orbit_mean(compute.orbit_mean(G,2),1);
        
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