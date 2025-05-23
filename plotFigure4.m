function plotFigure4(Parray, R, S, Narray)
%PLOTFIGURE4 This method plots figure 4 from the paper

%% Load in data

N = numel(Narray);
P = numel(Parray);
orbit_lsw_tikhonov_error= nan(P,R,N,S);
orbit_lsw_convex1_error = nan(P,R,N,S);
orbit_lsw_convex2_error = nan(P,R,N,S);
orbit_markov_error = nan(P,R,N,S);
orbit_pot_error    = nan(P,R);

str = '';
for s = 1:S 
    
    obj=load(sprintf('localdata/predictions/errors%g.mat',s));

    if s==1
        orbit_pot_error = max(obj.orbit_pot_error(:,:,1:end-3),[],3); % *_error(:,:,end-2:end) is the lyapunov exps and KY dim. err. Lets throw away this obs to compute E_max over \mathcal{B}
        orbit_pot_error = orbit_pot_error(:,1); % look only at the ordered library P_r = {1,...,P}. 
    end

    orbit_lsw_tikhonov_error(:,:,:,s)= max(obj.orbit_lsw_tikhonov_error(:,:,:,1:end-3),[],4); 
    orbit_lsw_convex1_error(:,:,:,s) = max(obj.orbit_lsw_convex1_error(:,:,:,1:end-3),[],4);
    orbit_lsw_convex2_error(:,:,:,s) = max(obj.orbit_lsw_convex2_error(:,:,:,1:end-3),[],4);

    orbit_markov_error(:,:,:,s) = max(obj.orbit_markov_error(:,:,:,1:end-3),[],4);

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
ybnds = [-4 0];
ylabels ={};
for i = ybnds(1):ybnds(2)
    ylabels{i-ybnds(1)+1} = ['$10^{',num2str(i),'}$'];
end

complete_libraries = [1 3 6 12 21 39 69 125];

%% define plotting parameters
clrs = lines(7);
palette = { '#DC267F','#785EF0','#648FFF','#FE6100','#FFB000'};
% palette = palette([3 5 2 1]);
ybnds = [-3.5 0];


% plot panel f (accuracy comparison)

% format figure
figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0]);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')

[~,n] = min(abs(Narray-10^6));

% plot_center_and_spread(Parray,orbit_pot_error,palette{1},.7,':');
plot_center_and_spread(Parray,permute(orbit_markov_error(:,:,n,:),[1 2 4 3]),palette{2},.7,'-');

plot_center_and_spread(Parray,permute(orbit_lsw_tikhonov_error(:,:,n,:),[1 2 4 3]),palette{3},.5,'-');
plot_center_and_spread(Parray,permute(orbit_lsw_convex1_error(:,:,n,:),[1 2 4 3]),palette{4},.5,'-');
plot_center_and_spread(Parray,permute(orbit_lsw_convex2_error(:,:,n,:),[1 2 4 3]),palette{5},.5,'-');

scatter(complete_libraries,orbit_pot_error(complete_libraries),110,'d','filled','CData',hex2rgb(palette{1}),'LineWidth',2);

% format axis

set(gca,'YScale','log','Xscale','log');
set(gcf,'color','w');
xlabel('$P$','Interpreter','latex');
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',24,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
ylabel('$E_\textrm{max}$','Interpreter','latex');
ylim(10.^ybnds);
yticks(10.^(-4:-1))
yticklabels(ylabels)
xticks(10.^(0:4))
grid on
box on
axis square

exportgraphics(gcf,'media/fig4f.pdf','ContentType','vector');

%% plot panels of weights at $P=125$ over each method

weights = cell(5,1);
p = 125;
r = 1;
n = 6;
obj=load('localdata/orbits/pot/weights.mat');
weights{1}=obj.w{p}(:,r);
obj=load('localdata/orbits/markov/weights1.mat');
weights{2}=obj.w{p}(:,r,n);
obj=load('localdata/orbits/lsw/weights1.mat');
weights(3:5) = {obj.w_tikhonov{p}(:,r,n), obj.w_convex1{p}(:,r,n), obj.w_convex2{p}(:,r,n)};
clear obj;
ytiks = [
    0 0.02 0.04;
    0 0.02 0.04;
    -4 0 4;
    0 0.3 0.6;
    0 0.02 0.04;
];
labels = 'deabc';

stability = zeros(p,1);
for i = 1:p
    orbit = load(sprintf('data/orbits/orbit%g.mat',i));
    stability(i) = orbit.floquetexponent;
end
[stability,ordering] = sort(stability,'ascend');

% panel a
for i = 1:numel(labels)

figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0],'color','w');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')

w = weights{i}(ordering);
l = stability(abs(w)>0);
w = w(abs(w)>0);
% w = log10(abs(w));
for j = 1:numel(w)
    plot(l(j)*[1 1],[0 w(j)],'color',hex2rgb(palette{i}),'LineWidth',2)
    hold on
end
scatter(l,w,100,'filled','CData',hex2rgb(palette{i}))
yline(0,'k-');
% b=bar(weights{i}(ordering));
% b.FaceColor = hex2rgb(palette{i});
% b.EdgeColor = hex2rgb(palette{i});

if i == 5
coeff = polyfit(l,w,1);
what = polyval(coeff,[0.78 1]);
plot([0.78 1],what,'k--','LineWidth',2);
text(.9, .025,['$w_p = ',sprintf('%.3f',coeff(1)),'\lambda_p^1 + ',sprintf('%.3f',coeff(2)),'$'],'interpreter','latex','FontSize',24,'HorizontalAlignment','center')
end

xlabel('$\lambda_p^1$');
% xticks([1 p]);
ylabel('$w_p$');
axis square
ylim(ytiks(i,[1 end]));
xlim([0.78 1])
yticks(ytiks(i,:))
set(gca,'fontsize',24,'YMinorGrid','on');


exportgraphics(gcf,sprintf('media/fig4%s.pdf',labels(i)),'ContentType','vector');

end

%% plot panels of weights over P over each method

weights = cell(5,1);
for i = 1:5
weights{i} = nan(125);
end
r = 1;
n = 6;
obj=load('localdata/orbits/pot/weights.mat');
complete_libraries = [1 3 6 12 21 39 69 125];
for p = 1:125
    if ismember(p,[1 3 6 12 21 39 69 125])
        weights{1}(p,1:p)=obj.w{p}(:,r);
        weights{1}(p,1:p)=weights{1}(p,1:p)/max(abs(weights{1}(p,1:p)));
    else
        q = find(p>complete_libraries,1,'last');
        q = complete_libraries(q);
        weights{1}(p,1:q)=weights{1}(q,1:q);
        weights{1}(p,q+1:p)=0;
    end
end
obj=load('localdata/orbits/markov/weights1.mat');
for p = 1:125
weights{2}(p,1:p)=obj.w{p}(:,r,n);
weights{2}(p,1:p)=weights{2}(p,1:p)/max(abs(weights{2}(p,1:p)));
end
obj=load('localdata/orbits/lsw/weights1.mat');
for p = 1:125
weights{3}(p,1:p)=obj.w_tikhonov{p}(:,r,n);
weights{4}(p,1:p)=obj.w_convex1{p}(:,r,n);
weights{5}(p,1:p)=obj.w_convex2{p}(:,r,n);
weights{3}(p,1:p)=weights{3}(p,1:p)/max(abs(weights{3}(p,1:p)));
weights{4}(p,1:p)=weights{4}(p,1:p)/max(abs(weights{4}(p,1:p)));
weights{5}(p,1:p)=weights{5}(p,1:p)/max(abs(weights{5}(p,1:p)));
end
clear obj;

labels = 'jkghi';
for i = 1:numel(labels)

figure
set(gcf,'defaultAxesColorOrder',[[0 0 0];0 0 0],'color','w');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')

% draw gray background denoting invalid cells
fill([0 P+1 P+1 0 0],[0 0 P+1 P+1 0],[1 1 1]*.8,'EdgeColor','none');
hold on
imagesc(weights{i}','AlphaData',~isnan(weights{i}'));
xticks([1 P]);
ylabel('$\hat{w}_p$');
yticks([1 P]);
xlabel('$P$');
clim([-1 1]);
colormap(rwb)
cb=colorbar;
cb.TickLabelInterpreter='latex';
set(gca,'Fontsize',24)
axis equal
xlim([0 P]+1/2);
ylim([0 P]+1/2);

exportgraphics(gcf,sprintf('media/fig4%s.pdf',labels(i)),'ContentType','vector');

end


end

function plot_center_and_spread(x,y,clr,alpha,ls)

y = reshape(y,size(y,1),[]);

% plot the center of y over samples
plot(x,median(y,2),'color',clr,'LineWidth',2,'LineStyle',ls);
hold on

%plot quantiles of y, if the number of chaotic samples is larger than one
if size(y,2)>1
    Q = quantile(y,[.25 .75],2);
    fill([x(:); flip(x(:))],[Q(:,1); flip(Q(:,2))],'k','Facecolor',clr,'FaceAlpha',alpha/5,'edgecolor',clr,'EdgeAlpha',alpha);
end

end