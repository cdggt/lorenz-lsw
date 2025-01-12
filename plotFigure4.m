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
        orbit_pot_error = max(obj.orbit_pot_error(:,:,1:end-1),[],3);
    end

    orbit_lsw_tikhonov_error(:,:,:,s)= max(obj.orbit_lsw_tikhonov_error(:,:,:,1:end-1),[],4);
    orbit_lsw_convex1_error(:,:,:,s) = max(obj.orbit_lsw_convex1_error(:,:,:,1:end-1),[],4);
    orbit_lsw_convex2_error(:,:,:,s) = max(obj.orbit_lsw_convex2_error(:,:,:,1:end-1),[],4);

    orbit_markov_error(:,:,:,s) = max(obj.orbit_markov_error(:,:,:,1:end-1),[],4);

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

plot_center_and_spread(Parray,orbit_pot_error,palette{1},.7,':');
plot_center_and_spread(Parray,permute(orbit_markov_error(:,:,n,:),[1 2 4 3]),palette{2},.7,'--');

plot_center_and_spread(Parray,permute(orbit_lsw_tikhonov_error(:,:,n,:),[1 2 4 3]),palette{3},.5,'-');
plot_center_and_spread(Parray,permute(orbit_lsw_convex1_error(:,:,n,:),[1 2 4 3]),palette{4},.5,'-');
plot_center_and_spread(Parray,permute(orbit_lsw_convex2_error(:,:,n,:),[1 2 4 3]),palette{5},.5,'-');

% format axis

set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
xlabel('$P$','Interpreter','latex');
set(gca,'YScale','log','Xscale','log');
set(gca,'Fontsize',20,'YMinorGrid','off','XMinorGrid','off');
set(gcf,'color','w');
ylabel('$E_\textrm{max}$','Interpreter','latex');
ylim(10.^ybnds);
yticks(10.^(-4:0))
xticks(10.^(0:4))
grid on
box on
axis square

exportgraphics(gcf,'media/fig4f.pdf','ContentType','vector');

%% plot panels of weights over each method

weights = cell(5,1);
p = 125;
s = 1;
n = 6;
obj=load('localdata/orbits/pot/weights.mat');
weights{1}=obj.w{p}(1:p);
obj=load('localdata/orbits/markov/weights1.mat');
weights{2}=obj.w{p}(:,s,n);
obj=load('localdata/orbits/lsw/weights1.mat');
weights(3:5) = {obj.w_tikhonov{p}(:,s,n), obj.w_convex1{p}(:,s,n), obj.w_convex2{p}(:,s,n)};
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
text(.9, .025,['$w_p = ',sprintf('%.3f',coeff(1)),'\lambda_p + ',sprintf('%.3f',coeff(2)),'$'],'interpreter','latex','FontSize',24,'HorizontalAlignment','center')
end

xlabel('$\lambda_p$');
% xticks([1 p]);
ylabel('$w_p$');
axis square
ylim(ytiks(i,[1 end]));
xlim([0.78 1])
yticks(ytiks(i,:))
set(gca,'fontsize',24,'YMinorGrid','on');


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