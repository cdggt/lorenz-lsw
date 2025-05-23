function table = plotTable1(Parray, R, S, Narray)

%% load in data of table 

[~,pindex] = min(abs(Parray-21));
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

        orbit_uniform_error   = squeeze(obj.orbit_uniform_error(pindex,:,:)); % P x R x Nobs
        orbit_pot_error       = squeeze(obj.orbit_pot_error(pindex,1,:));  % P x R x Nobs
        snippet_uniform_error = squeeze(obj.snippet_uniform_error(pindex,:,:));  % P x R x Nobs

        table(:,1) = orbit_pot_error;
        table(:,2) = median(orbit_uniform_error,1)';
        table(:,5) = median(snippet_uniform_error,1)';

    end

    orbit_markov_error(:,:,s) = squeeze(obj.orbit_markov_error(pindex,:,nindex,:));
    snippet_markov_error(:,:,s) = squeeze(obj.snippet_markov_error(pindex,:,nindex,:));

    orbit_lsw_error(:,:,s)    = squeeze(obj.orbit_lsw_tikhonov_error(pindex,:,nindex,:));
    snippet_lsw_error(:,:,s)    = squeeze(obj.snippet_lsw_tikhonov_error(pindex,:,nindex,:));
    
    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',s,S);
    fprintf(str);

end

table(:,3) = median(reshape(permute(orbit_markov_error,[2 1 3]),Nobs,[]),2);
table(:,4) = median(reshape(permute(orbit_lsw_error,[2 1 3]),Nobs,[]),2);
table(:,6) = median(reshape(permute(snippet_markov_error,[2 1 3]),Nobs,[]),2);
table(:,7) = median(reshape(permute(snippet_lsw_error,[2 1 3]),Nobs,[]),2);

Erel = log10(table);

%% plot table

figure;
setlatexlabels

heatmap(round(Erel,1));
ax = gca;
ax.YData = {"$1$","$x$","$y$","$z$","$x^2$","$xy$","$xz$","$y^2$","$yz$","$z^2$","$\lambda^1$","$\lambda^3$","$d_{KY}$"};
ax.XData = {"POT (O)", "Uniform (O)","Markov (O)","LSW (O)","Markov (S)","Uniform (S)","LSW (S)"};
ax.Title = '$\log(E_\textrm{rel})$';
ax.Interpreter='latex';
colormap(summer)

set(gcf,'color','w');
set(gca,'fontsize',14)

exportgraphics(gcf,'media/tab1.pdf','ContentType','vector');

%% print table 

labels = {"$1$","$x$","$y$","$z$","$x^2$","$xy$","$xz$","$y^2$","$yz$","$z^2$","$\lambda^1$","$\lambda^3$","$d_{KY}$"};
for i = 1:numel(labels)
    fprintf('%s & %s & %s & %s & %s & %s & %s & %s \\\\ \n',labels{i},format(Erel(i,1)),format(Erel(i,2)),format(Erel(i,3)),format(Erel(i,4)),format(Erel(i,5)),format(Erel(i,6)),format(Erel(i,7)));
end

end

function x = format(x)
    if isnumeric(x)

        if isnan(x)
            x = '$-$';
        else
            if x<-12
                x = '$-\boldsymbol{\infty}$';
            else
                x = ['$',sprintf('%.1f',round(x,1)),'$'];
            end
        end
    else
        error('type not recognized');
    end
    x = ['\new{',x,'}'];
end