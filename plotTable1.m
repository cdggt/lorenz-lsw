function table = plotTable1(Parray, R, S, Narray)

%% load in data of table 

Nobs = 11;
[~,pindex] = min(abs(Parray-25));
[~,nindex] = min(abs(Narray-10^6));

str = '';
for s = 1:S 
    
    obj=load(sprintf('localdata/predictions/errors%g.mat',s));

    if s==1

        Nobs = size(obj.orbit_markov_error,4);
        orbit_markov_error = nan(R,Nobs,S);
        snippet_markov_error = nan(R,Nobs,S);
        orbit_lsw_error    = nan(R,Nobs,S);
        snippet_lsw_error    = nan(R,Nobs,S);
        table = zeros(Nobs,7);

        % compute errors for given p. Results are [1 x R x Nobs]
        orbit_uniform_error   = obj.orbit_uniform_error(pindex,:,:);
        orbit_pot_error       = obj.orbit_pot_error(pindex,:,:);
        snippet_uniform_error = obj.snippet_uniform_error(pindex,:,:);

        % make results [Nobs x R]
        orbit_uniform_error   = reshape(permute(orbit_uniform_error,[3 2 1]),Nobs,[]);
        orbit_pot_error       = reshape(permute(orbit_pot_error,[3 2 1]),Nobs,[]);
        snippet_uniform_error = reshape(permute(snippet_uniform_error,[3 2 1]),Nobs,[]);

        % take median over R
        table(:,1) = median(orbit_pot_error,2);
        table(:,2) = median(orbit_uniform_error,2);
        table(:,5) = median(snippet_uniform_error,2);

    end

    orbit_markov_error(:,:,s)   = permute(obj.orbit_markov_error(pindex,:,nindex,:),[2 4 1 3]);
    snippet_markov_error(:,:,s) = permute(obj.snippet_markov_error(pindex,:,nindex,:),[2 4 1 3]);

    orbit_lsw_error(:,:,s)      = permute(obj.orbit_lsw_tikhonov_error(pindex,:,nindex,:),[2 4 1 3]);
    snippet_lsw_error(:,:,s)    = permute(obj.snippet_lsw_tikhonov_error(pindex,:,nindex,:),[2 4 1 3]);
    
    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',s,S);
    fprintf(str);

end

table(:,3) = median(reshape(orbit_markov_error,Nobs,[]),2);
table(:,4) = median(reshape(orbit_lsw_error,Nobs,[]),2);
table(:,6) = median(reshape(snippet_markov_error,Nobs,[]),2);
table(:,7) = median(reshape(snippet_lsw_error,Nobs,[]),2);

Erel = log10(table);

%% plot table

figure;
setlatexlabels

heatmap(round(Erel,1));
ax = gca;
ax.YData = {"$1$","$x$","$y$","$z$","$x^2$","$xy$","$xz$","$y^2$","$yz$","$z^2$","$\lambda$"};
ax.XData = {"POT${}_{orbits}$", "Uniform${}_{orbits}$","Markov${}_{orbits}$","LSW${}_{orbits}$","Markov${}_{snippets}$","Uniform${}_{snippets}$","LSW${}_{snippets}$"};
ax.Title = '$\log(E_\textrm{rel})$';
ax.Interpreter='latex';
colormap(summer)

set(gcf,'color','w');
set(gca,'fontsize',14)

exportgraphics(gcf,'media/tab1.pdf','ContentType','vector');


end