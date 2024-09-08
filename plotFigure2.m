function plotFigure2(Parray, R, S, Narray)
%PLOTFIGURE2 This method plots figure 2 from the paper

%% Load in data

N = numel(Narray);
P = numel(Parray);
orbit_lsw_error    = nan(P,R,N,S);
orbit_markov_error = nan(P,R,N,S);
orbit_uniform_error= nan(P,R);
orbit_pot_error    = nan(P,R);
snippet_lsw_error    = nan(P,R,N,S);
snippet_markov_error = nan(P,R,N,S);
snippet_uniform_error= nan(P,R);
str = '';
for s = 1:S 
    
    obj=load(sprintf('localdata/predictions/errors%g.mat',s));

    if s==1
        orbit_uniform_error   = max(obj.orbit_uniform_error(:,:,1:end-1),[],3);
        orbit_pot_error       = max(obj.orbit_pot_error(:,:,1:end-1),[],3);
        snippet_uniform_error = max(obj.snippet_uniform_error(:,:,1:end-1),[],3);
    end

    orbit_lsw_error(:,:,:,s)    = max(obj.orbit_lsw_error(:,:,:,1:end-1),[],4);
    orbit_markov_error(:,:,:,s) = max(obj.orbit_markov_error(:,:,:,1:end-1),[],4);

    snippet_lsw_error(:,:,:,s)    = max(obj.snippet_lsw_error(:,:,:,1:end-1),[],4);
    snippet_markov_error(:,:,:,s) = max(obj.snippet_markov_error(:,:,:,1:end-1),[],4);

    fprintf(repmat('\b',1,numel(str)));
    str = sprintf('\t %g / %g \n',s,S);
    fprintf(str);

end

%% plot 

line_xcrds = {
    Narray(1:6),
    Narray(1:6),
    Narray(1:6),
    Narray(1:6)
};
line_ycrds = {
    10^(-0.5)*line_xcrds{1}.^(-1/2),
    10^(0.25)*line_xcrds{2}.^(-1/2),
    10^(-0.15)*line_xcrds{3}.^(-1/2),
    10^(0.40)*line_xcrds{4}.^(-1/2)
};

%% define plotting parameters

% palette = {cmap(2,:),cmap(3,:),'#1b9e77',[0 0 0],cmap(4,:)};
% palette = {'#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'};
palette = {'#648FFF','#785EF0','#DC267F','#FE6100','#FFB000'};
% palette = {'#0077BB','#009988','#000000','#CC3311','#EE3377'};
% palette = {cmap(1,:),cmap(2,:),'#000000',cmap(4,:),cmap(5,:)};
palette = palette([3 5 2 1]);
ybnds = [-4 0];
ylabels ={};
for i = ybnds(1):ybnds(2)
    ylabels{i-ybnds(1)+1} = ['$10^{',num2str(i),'}$'];
end

%% plot panel a

% format figure
figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0]);
pos = get(gcf,'Position');
pos = set(gcf,'Position',[pos(1) pos(2) pos(3)*.9 pos(4)*.7]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')

% plot performance over N for fixed P
[~,p] = min(abs(Parray-5));
subplot(1,2,1);

plot_center_and_spread(Narray,repmat(orbit_pot_error(p,:),[numel(Narray),1]),palette{1},.7,':');
plot_center_and_spread(Narray,repmat(orbit_uniform_error(p,:),[numel(Narray),1]),palette{2},.7,'-.');
plot_center_and_spread(Narray,permute(orbit_markov_error(p,:,:,:),[3 2 4 1]),palette{3},.7,'--');
plot_center_and_spread(Narray,permute(orbit_lsw_error(p,:,:,:),[3 2 4 1]),palette{4},.5,'-');

plot(line_xcrds{1},line_ycrds{1},'k-','LineWidth',2)

% format axis
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.1700    0.200    0.35    0.75])
ylabel('$E_\textrm{max}$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylim(10.^ybnds);
xlim(10.^[1 6])
xticks(10.^(1:2:6))
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels(ylabels)
grid on
box on
xcrd = 10^(log10(Narray(1))+.05*diff(log10(Narray([1 end]))));
text(xcrd,10^(-3.6),['$P=',num2str(Parray(p)),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

% plot performance over N for fixed P
[~,p] = min(abs(Parray-125));
subplot(1,2,2);

plot_center_and_spread(Narray,repmat(orbit_pot_error(p,:),[numel(Narray),1]),palette{1},.7,':');
plot_center_and_spread(Narray,repmat(orbit_uniform_error(p,:),[numel(Narray),1]),palette{2},.7,'-.');
plot_center_and_spread(Narray,permute(orbit_markov_error(p,:,:,:),[3 2 4 1]),palette{3},.7,'--');
plot_center_and_spread(Narray,permute(orbit_lsw_error(p,:,:,:),[3 2 4 1]),palette{4},.5,'-');

plot(line_xcrds{2},line_ycrds{2},'k-','LineWidth',2)

% format axis
% ylabel('$E_\textrm{max}$','Interpreter','latex');
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.5703    0.200    0.35    0.75])
xlabel('$N$','Interpreter','latex');
ylim(10.^ybnds);
xlim(10.^[1 6])
xticks(10.^(1:2:6))
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels({'','','','',''});
grid on
box on
text(xcrd,10^(-3.6),['$P=',num2str(Parray(p)),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

exportgraphics(gcf,'media/fig2a.pdf','ContentType','vector');

%% plot panel b

% format figure
figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0]);
pos = get(gcf,'Position');
pos = set(gcf,'Position',[pos(1) pos(2) pos(3)*.9 pos(4)*.7]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
cmap = lines(8);

% plot performance over P for fixed N
[~,n] = min(abs(Narray-10^3));
subplot(1,2,1);

plot_center_and_spread(Parray,orbit_pot_error,palette{1},.7,':');
plot_center_and_spread(Parray,orbit_uniform_error,palette{2},.7,'-.');
plot_center_and_spread(Parray,permute(orbit_markov_error(:,:,n,:),[1 2 4 3]),palette{3},.7,'--');
plot_center_and_spread(Parray,permute(orbit_lsw_error(:,:,n,:),[1 2 4 3]),palette{4},.5,'-');

% format axis
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.1700    0.200    0.35    0.75])
ylabel('$E_\textrm{max}$','Interpreter','latex');
xlabel('$P$','Interpreter','latex');
ylim(10.^ybnds);
xticks(10.^(0:4))
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels(ylabels)
grid on
box on
xcrd = 10^(log10(Parray(1))+.05*diff(log10(Parray([1 end]))));
text(xcrd,10^(-3.6),['$N=10^',num2str(log10(Narray(n))),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

% plot performance over P for fixed N
[~,n] = min(abs(Narray-10^6));
subplot(1,2,2);

plot_center_and_spread(Parray,orbit_pot_error,palette{1},.7,':');
plot_center_and_spread(Parray,orbit_uniform_error,palette{2},.7,'-.');
plot_center_and_spread(Parray,permute(orbit_markov_error(:,:,n,:),[1 2 4 3]),palette{3},.7,'--');
plot_center_and_spread(Parray,permute(orbit_lsw_error(:,:,n,:),[1 2 4 3]),palette{4},.5,'-');

% format axis
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.5703    0.200    0.35    0.75])
xlabel('$P$','Interpreter','latex');
ylim(10.^ybnds);
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels({'','','','',''});
xticks(10.^(0:4))
grid on
box on
text(xcrd,10^(-3.6),['$N=10^',num2str(log10(Narray(n))),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

exportgraphics(gcf,'media/fig2b.pdf','ContentType','vector');

%% plot panel c

% format figure
figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0]);
pos = get(gcf,'Position');
pos = set(gcf,'Position',[pos(1) pos(2) pos(3)*.9 pos(4)*.7]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')

% plot performance over N for fixed P
[~,p] = min(abs(Parray-5));
subplot(1,2,1);

plot_center_and_spread(Narray,repmat(snippet_uniform_error(p,:),[numel(Narray),1]),palette{2},.7,'-.');
plot_center_and_spread(Narray,permute(snippet_markov_error(p,:,:,:),[3 2 4 1]),palette{3},.7,'--');
plot_center_and_spread(Narray,permute(snippet_lsw_error(p,:,:,:),[3 2 4 1]),palette{4},.5,'-');

plot(line_xcrds{3},line_ycrds{3},'k-','LineWidth',2)

% format axis
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.1700    0.200    0.35    0.75])
ylabel('$E_\textrm{max}$','Interpreter','latex');
xlabel('$N$','Interpreter','latex');
ylim(10.^ybnds);
xlim(10.^[1 6])
xticks(10.^(1:2:6))
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels(ylabels)
grid on
box on
xcrd = 10^(log10(Narray(1))+.05*diff(log10(Narray([1 end]))));
text(xcrd,10^(-3.6),['$P=',num2str(Parray(p)),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

% plot performance over N for fixed P
[~,p] = min(abs(Parray-125));
subplot(1,2,2);

plot_center_and_spread(Narray,repmat(snippet_uniform_error(p,:),[numel(Narray),1]),palette{2},.7,'-.');
plot_center_and_spread(Narray,permute(snippet_markov_error(p,:,:,:),[3 2 4 1]),palette{3},.7,'--');
plot_center_and_spread(Narray,permute(snippet_lsw_error(p,:,:,:),[3 2 4 1]),palette{4},.5,'-');

plot(line_xcrds{4},line_ycrds{4},'k-','LineWidth',2)

% format axis
% ylabel('$E_\textrm{max}$','Interpreter','latex');
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.5703    0.200    0.35    0.75])
xlabel('$N$','Interpreter','latex');
ylim(10.^ybnds);
xlim(10.^[1 6])
xticks(10.^(1:2:6))
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels({'','','','',''});
grid on
box on
text(xcrd,10^(-3.6),['$P=',num2str(Parray(p)),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

exportgraphics(gcf,'media/fig2c.pdf','ContentType','vector');

%% plot panel d

% format figure
figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0]);
pos = get(gcf,'Position');
pos = set(gcf,'Position',[pos(1) pos(2) pos(3)*.9 pos(4)*.7]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')
cmap = lines(8);

% plot performance over P for fixed N
[~,n] = min(abs(Narray-10^3));
subplot(1,2,1);

plot_center_and_spread(Parray,snippet_uniform_error,palette{2},.7,'-.');
plot_center_and_spread(Parray,permute(snippet_markov_error(:,:,n,:),[1 2 4 3]),palette{3},.7,'--');
plot_center_and_spread(Parray,permute(snippet_lsw_error(:,:,n,:),[1 2 4 3]),palette{4},.5,'-');

% format axis
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.1700    0.200    0.35    0.75])
ylabel('$E_\textrm{max}$','Interpreter','latex');
xlabel('$P$','Interpreter','latex');
ylim(10.^ybnds);
xticks(10.^(0:4))
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels(ylabels)
grid on
box on
xcrd = 10^(log10(Parray(1))+.05*diff(log10(Parray([1 end]))));
text(xcrd,10^(-3.6),['$N=10^',num2str(log10(Narray(n))),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

% plot performance over P for fixed N
[~,n] = min(abs(Narray-10^6));
subplot(1,2,2);

plot_center_and_spread(Parray,snippet_uniform_error,palette{2},.7,'-.');
plot_center_and_spread(Parray,permute(snippet_markov_error(:,:,n,:),[1 2 4 3]),palette{3},.7,'--');
plot_center_and_spread(Parray,permute(snippet_lsw_error(:,:,n,:),[1 2 4 3]),palette{4},.5,'-');

% format axis
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
set(gca,'Position',[0.5703    0.200    0.35    0.75])
xlabel('$P$','Interpreter','latex');
ylim(10.^ybnds);
yticks(10.^(ybnds(1):ybnds(2)))
yticklabels({'','','','',''});
xticks(10.^(0:4))
grid on
box on
text(xcrd,10^(-3.6),['$N=10^',num2str(log10(Narray(n))),'$'],'Interpreter','latex','FontSize',20,'BackgroundColor','w','EdgeColor','k','HorizontalAlignment','left');

exportgraphics(gcf,'media/fig2d.pdf','ContentType','vector');

end

function plot_center_and_spread(x,y,clr,alpha,ls)

y = reshape(y,size(y,1),[]);

% plot the center of y over samples
plot(x,mean(y,2),'color',clr,'LineWidth',2,'LineStyle',ls);
hold on

% plot quantiles of y, if the number of chaotic samples is larger than one
if size(y,2)>1
    Q = quantile(y,[.25 .75],2);
    fill([x(:); flip(x(:))],[Q(:,1); flip(Q(:,2))],'k','Facecolor',clr,'FaceAlpha',alpha/5,'edgecolor',clr,'EdgeAlpha',alpha);
end

end