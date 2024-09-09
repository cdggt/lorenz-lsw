function plotFigure1(recompute)

filename = './localdata/fig1.mat';
if ~isfile(filename)||recompute

    theta = logspace(-1,5,7*4); % the range of theta to investigate
    N = numel(theta);
    P = 30; % the library size to investiage

    Ko = zeros(P,P,N); % the orbit correlation matrix, K_pq^o
    Ks = zeros(P,P,N); % the snippet correlation matrix, K_pq^s
    K = zeros(P); % a temporary matrix to let the PARFOR loop to work

    deyeo = zeros(N,1); % the distance of K_pq^o to the identity matrix
    deyes = zeros(N,1); % the distance of K_pq^s to the identity matrix
    doneo = zeros(N,1); % the distance of K_pq^o to the matrix of ones
    dones = zeros(N,1); % the distance of K_pq^s to the matrix of ones

    On = ones(P); % the identity matrix
    Id = eye(P); % the matrix of ones


    %% compute data

    fprintf('computing data for figure 1...\n');

    str='';
    for i=1:N

        th = theta(i);

        % orbits
        for p = 1:P

            % load in data of orbit p
            orbitp =load(sprintf('./data/orbits/orbit%g.mat',p));
            parfor q = p:P

                % load in data of orbit p
                orbitq =load(sprintf('./data/orbits/orbit%g.mat',q));

                % compute the integral of the Gaussian kernel
                distance = (orbitp.x-orbitq.x').^2+(orbitp.y-orbitq.y').^2+(orbitp.z-orbitq.z').^2;
                G = exp(-1/(4*th).*distance);
                K(p,q) = compute.orbit_mean(compute.orbit_mean(G,2),1);

            end

        end
        % compute Frobenius distances
        Ko(:,:,i) = K+K'-diag(diag(K));
        deyeo(i) = norm(Ko(:,:,i)-Id);
        doneo(i) = norm(Ko(:,:,i)-On);

        % snippets
        for p = 1:P

            % load in data of snippet p
            snippetp =load(sprintf('./localdata/snippets/snippet%g.mat',p));
            parfor q = p:P

                % load in data of snippet p
                snippetq =load(sprintf('./localdata/snippets/snippet%g.mat',q));

                % compute the integral of the Gaussian kernel
                distance = (snippetp.x-snippetq.x').^2+(snippetp.y-snippetq.y').^2+(snippetp.z-snippetq.z').^2;
                G = exp(-1/(4*th).*distance);
                K(p,q) = compute.snippet_mean(compute.snippet_mean(G,2),1);

            end

        end
        % compute Frobenius distances
        Ks(:,:,i) = K+K'-diag(diag(K));
        deyes(i) = norm(Ks(:,:,i)-Id);
        dones(i) = norm(Ks(:,:,i)-On);

        fprintf(repmat('\b',1,numel(str)));
        str = sprintf('\t %g / %g \n',i,N);
        fprintf(str);

        save(filename,'theta','P','doneo','deyeo','dones','deyes');

    end

else

    load(filename,'theta','P','doneo','deyeo','dones','deyes');

end

%% plot data
figure
setlatexlabels

clr = {'#648FFF','#785EF0','#DC267F','#FE6100','#FFB000'};

p1 = plot(log10(theta),doneo/P,':','LineWidth',3,'color',clr{1});
hold on;
p2 = plot(log10(theta),deyeo/P,'LineWidth',2,'color',clr{4});
plot(log10(theta),dones/P,'o','Color',p1.Color,'LineWidth',2,'MarkerSize',10);
plot(log10(theta),deyes/P,'x','Color',p2.Color,'LineWidth',3,'MarkerSize',10);

xline(2,'k-','LineWidth',1);
yticks([0 1]);
yticklabels({'0','1'});
ylim([0 1]);
xlim([-1 5]);
xlabel('$\log_{10}\theta$','Interpreter','latex');
ylabel('$P^{-1}\| \ \cdot \ \|_\textrm{F}$','Interpreter','latex');
set(gcf,'color','w');
set(gca,'fontsize',30);

exportgraphics(gcf,'media/fig1.pdf');

end


