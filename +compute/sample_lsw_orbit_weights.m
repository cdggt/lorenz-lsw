function sample_lsw_orbit_weights(recompute,sampleNumber,Parray,Narray,theta)
%SAMPLE_LSW_ORBIT_WEIGHTS this method computes the LSW weights of orbits,
%for each P, R, and N, at a specific sample index S=sampleIndex. 

filename = sprintf('localdata/orbits/lsw/weights%g.mat',sampleNumber);
if isfile(filename)&&~recompute

    fprintf('orbit lsw weights for sample %g already exists. \n',sampleNumber);

else

    fprintf('computing orbit lsw averages from the sample number %g...\n',sampleNumber)

    load('data/library_permutations.mat','permutations');

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
        orbit =load(sprintf('data/orbits/orbit%g.mat',p));

        D = pdist2([orbit.x orbit.y orbit.z], sample);
        ap = compute.orbit_mean(exp(-1/(4*theta).*D.^2),1);

        ap = cumsum(ap);
        averages(p,:) = ap(Narray)./Narray;

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',p,P);
        fprintf(str);

    end

    fprintf('computing orbit lsw weights from the sample number %g...\n',sampleNumber)

    load('localdata/orbits/lsw/correlations.mat','K');
    R = size(permutations,2);
    N = numel(Narray);
    w_tikhonov = cell(P,1); % weights using Tihkonov regularization
    w_convex1 = cell(P,1); % weights using lsqnonneg
    w_convex2 = cell(P,1); % weights using fmincon
    ridgeparameter = 1e-6;
    str = '';
    for i = 1:P

        p = Parray(i);
        I = eye(p);
        w_tikhonov{i} = zeros(p,R,N);
        w_convex1{i} = zeros(p,R,N);
        w_convex2{i} = zeros(p,R,N);

        for r = 1:R

            ind = permutations(1:p,r);
            kernel = K(ind,ind);

            % tikhonov regularized weights
            w_tikhonov{i}(:,r,:) = (kernel+ridgeparameter*I) \ averages(ind,:);

            for n = 1:N

                % convex method 1
                w = lsqnonneg(kernel,averages(ind,n));
                w_convex1{i}(:,r,n) = w/sum(w);

                % convex method 2
                w_convex2{i}(:,r,n) = fminconvex(kernel,averages(ind,n));
         
            end
            
        end

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',i,P);
        fprintf(str);

    end

    % save out data
    save(filename,'w_tikhonov','w_convex1','w_convex2','theta');
    fprintf('saved results to `%s`\n',filename)

end

end

function w = fminconvex(A,b)

% Make initial guess
n = size(A, 2);
guess = ones(n,1)/n;

% Objective function including the l1 norm penalty
objective = @(x) 0.5 * norm(A * x - b)^2;

% Constraints: sum(w) = 1 and w >= 0
Aeq = ones(1, n);
beq = 1;
lb = zeros(n, 1);

% Solve using fmincon
options = optimoptions('fmincon', 'Algorithm', 'sqp','MaxFunctionEvaluations',1e3,'OptimalityTolerance',1e-6,'Display','off');
w = fmincon(objective, guess, [], [], Aeq, beq, lb, [], [], options);

end