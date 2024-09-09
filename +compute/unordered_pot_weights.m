function w = unordered_pot_weights(Pr)
%UNORDERED_POT_WEIGHTS this method computes pot weights of the orbits in
%P_r. See the supplemental material for description of what this code is
%computing.

if isscalar(Pr)
    w=1;
    return;
end

maxtopologicalperiod = 0;
for i = 1:numel(Pr)
    orbit = load(sprintf('data/orbits/orbit%g.mat',Pr(i)));
    maxtopologicalperiod = max(maxtopologicalperiod,orbit.topologicalperiod);
end

% run a newton solver to find s such that F(s,0) = 0, see ( Eq 4 of supp. material )
s = 0; % initial guess
maxIter = 1e5;
for i = 1:maxIter

    % evaluate F and its derivatives at this value of s
    [Q,dQ_ds,d2Q_dbeta_dmu] = compute_determinant_derivatives(Pr,s,maxtopologicalperiod);

    % compute newton step
    F = 1-sum(Q);
    dF_ds = sum(-dQ_ds);
    update = - F / dF_ds;

    % if update size is beneath tolerance, end iterations
    if abs(update) <  1e-8
        break
    end

    % otherwise, step
    s = s + update;

end

% print failure to reach optimum if necessary
if i == maxIter

    fprintf("POT failed to converge exactly!\n")

end

% compute weights using Eq 2 of supp. material
w = -sum(d2Q_dbeta_dmu,1)' / sum(dQ_ds);

end

function [Q,dQ_ds,d2Q_dbeta_dmu] = compute_determinant_derivatives(orbitIndices,s,maxtopologicalperiod)

P = numel(orbitIndices);

% allocate memory for arrays
C = zeros(1,maxtopologicalperiod);
dC_ds = zeros(1,maxtopologicalperiod);
d2C_dbeta_dmu = zeros(maxtopologicalperiod,P);
Q = zeros(1,maxtopologicalperiod); %
dQ_ds = zeros(1,maxtopologicalperiod);
d2Q_dbeta_dmu = zeros(maxtopologicalperiod,length(orbitIndices));

%% Compute trace coefficients C_j

% for each orbit p = 1,...,P
for p = 1:P

    orbit = load(sprintf('data/orbits/orbit%g.mat',orbitIndices(p)));
    period = orbit.period;

    % for each repeat, r, of orbit p, such that r*n_p = i (see Eq. 3 of supp. material)
    for r = 1:floor(maxtopologicalperiod/orbit.topologicalperiod)

        % compute stability of each orbit det(1-M_p^r). see just below Eq. 3 of supp. material
        stability = abs(prod(1-exp(period * orbit.floquetexponent)^r));
        j = r*orbit.topologicalperiod;
        weight = orbit.topologicalperiod * exp(-r * period * (s))/stability;

        % add term to trace coefficient sum (Eq 3 of supp. material)
        C(j) = C(j) + weight;

        % add term to trace coefficient derivative sum (Eq 6 of supp. material)
        dC_ds(j) = dC_ds(j) + weight * (-r * period);

        % add term to trace coefficient derivative sum (Eq 8 of supp. material)
        d2C_dbeta_dmu(j,p) = d2C_dbeta_dmu(j,p) + weight * period * (+r);

    end
end

%% Compute Cummulants Q_j

for k = 1:maxtopologicalperiod

    j = 1:(k-1);

    % compute Q_k from Eq 5 in supp. material
    Q(k) = (C(k) - sum(C(flip(j)) .* Q(j)))/k;

    % compute partial derivatives of Q_k via application of the chain rule on Eqs 5
    % with Eqs 6-8 in supp. material
    dQ_ds(k) = (dC_ds(k) - sum(dC_ds(flip(j)) .* Q(j)) - sum(C(flip(j)) .* dQ_ds(j)))/k;
    d2Q_dbeta_dmu(k,:) = (d2C_dbeta_dmu(k,:) - sum(d2C_dbeta_dmu(flip(j),:) .* Q(j)') - sum(C(flip(j))' .* d2Q_dbeta_dmu(j,:)))/k;

end

end