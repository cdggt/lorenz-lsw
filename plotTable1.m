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

        orbit_uniform_error   = squeeze(obj.orbit_uniform_error(pindex,:,:));
        orbit_pot_error       = squeeze(obj.orbit_pot_error(pindex,:,:));
        snippet_uniform_error = squeeze(obj.snippet_uniform_error(pindex,:,:));

        table(:,1) = median(orbit_pot_error,1)';
        table(:,2) = median(orbit_uniform_error,1)';
        table(:,3) = median(snippet_uniform_error,1)';

    end

    orbit_markov_error(:,:,s) = squeeze(obj.orbit_markov_error(pindex,:,nindex,:));
    snippet_markov_error(:,:,s) = squeeze(obj.snippet_markov_error(pindex,:,nindex,:));

    orbit_lsw_error(:,:,s)    = squeeze(obj.orbit_lsw_error(pindex,:,nindex,:));
    snippet_lsw_error(:,:,s)    = squeeze(obj.snippet_lsw_error(pindex,:,nindex,:));
    
    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',s,S);
    fprintf(str);

end

table(:,4) = median(reshape(permute(orbit_markov_error,[2 1 3]),Nobs,[]),2);
table(:,5) = median(reshape(permute(snippet_markov_error,[2 1 3]),Nobs,[]),2);
table(:,6) = median(reshape(permute(orbit_lsw_error,[2 1 3]),Nobs,[]),2);
table(:,7) = median(reshape(permute(snippet_lsw_error,[2 1 3]),Nobs,[]),2);

Erel = log10(table);

%% plot table

figure;
setlatexlabels

heatmap(round(Erel,2));
ax = gca;
ax.YData = {"$1$","$x$","$y$","$z$","$x^2$","$xy$","$xz$","$y^2$","$yz$","$z^2$","$\lambda$"};
ax.XData = {"POT$_{o}$", "Uniform$_{o}$","Uniform$_{s}$","Markov$_{o}$","Markov$_{s}$","LSW$_{o}$","LSW$_{s}$"};
ax.Title = '$\log(E_\textrm{rel})$';
ax.Interpreter='latex';
colormap(summer)

set(gcf,'color','w');
set(gca,'fontsize',14)

exportgraphics(gcf,'media/tab1.pdf','ContentType','vector');


end