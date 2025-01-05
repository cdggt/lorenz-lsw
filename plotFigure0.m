function plotFigure0(recompute)

filename = './localdata/fig0.mat';
if ~isfile(filename)||(nargin>0&&recompute)

    fprintf('computing data for figure 0...\n');

    % kernel variance 
    theta = 1^2; % produces a more visual/didactic plot than the optimal value, theta=10^2. 

    % decide grid over which to plot functions. Larger/Denser grids are more
    % computationally expensive. We are marginalizing over y, so it does
    % not appear in this script.
    x = linspace(-30, 30,600/2); 
    z = linspace(0, 50,500/2);
    [X,Z] = ndgrid(x,z);
    
    % Compute orbit density
    p = 19;
    orbit = load(sprintf('data/orbits/orbit%g.mat',p));
    orbit_trajectory = [orbit.x orbit.y orbit.z];
    T = size(orbit_trajectory,1);
    rho_p = zeros(numel(x),numel(z));
    for t = 1:T
        G = exp(-((orbit_trajectory(t,1)-X).^2+(orbit_trajectory(t,3)-Z).^2)/(4*theta));
        rho_p = rho_p+G;
    end

    % Compute chaotic density and histogram
    rho = zeros(numel(x),numel(z));
    histogram = zeros(numel(x),numel(z));
    for n = 1:50
        chaos = load(sprintf('localdata/chaos/sample%g.mat',n));
        chaotic_trajectory = cat(1,[chaos.x' chaos.y' chaos.z']);
        T = size(chaotic_trajectory,1);
        for t = 1:T

            G = exp(-((chaotic_trajectory(t,1)-X).^2+(chaotic_trajectory(t,3)-Z).^2)/(4*theta));
            rho = rho+G;

            dist = pdist2([X(:) Z(:)], [chaotic_trajectory(t,1) chaotic_trajectory(t,3)]);
            [~,i] = min(dist);
            histogram(i) = histogram(i)+1;

        end

    end
    

    % compute a sample chaotic trajectory to plot over the histogram

    % Compute chaotic density and histogram
    chaotic_trajectory = [5*rand; 5*rand; 20];
    timestep = 2e-3;
    for i = 1:1000
        chaotic_trajectory(:,1) = lorenz_rk4(chaotic_trajectory(:,1),timestep);
    end

    T = 2267; % this makes the chaotic snippet the same length as orbit 19
    for t = 2:T
        chaotic_trajectory(:,t) = lorenz_rk4(chaotic_trajectory(:,t-1),timestep);
    end
    chaotic_trajectory = chaotic_trajectory';

    save(filename,'x','z','rho','rho_p','histogram','orbit_trajectory','chaotic_trajectory');

else

    load(filename,'x','z','rho','rho_p','histogram','orbit_trajectory','chaotic_trajectory');

end

%% Compute chaotic density and histogram
chaotic_trajectory = [5*rand; 5*rand; 20];
timestep = 2e-3;
for i = 1:1000
    chaotic_trajectory(:,1) = lorenz_rk4(chaotic_trajectory(:,1),timestep);
end

T = 1500; % this makes the chaotic snippet the same length as orbit 19
for t = 2:T
    chaotic_trajectory(:,t) = lorenz_rk4(chaotic_trajectory(:,t-1),timestep);
end
chaotic_trajectory = chaotic_trajectory';

%% plot panels

lw = 2;
fs = 24;

figure
setlatexlabels
imagesc(rho_p',XData=x,YData=z); 
set(gca,'ydir','normal');
colormap(flip(bone));
hold on
plot(orbit_trajectory([1:end 1],1),orbit_trajectory([1:end 1],3),'r','LineWidth',lw);
xlabel('$x$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
set(gcf,'color','w');
set(gca,'fontsize',fs);
axis equal
ylim([0 50]);
xlim([-25 25]);
xticks([-20 0 20]);
yticks([0 25 50]);

exportgraphics(gcf,'media/fig0a.pdf');

dir = diff(chaotic_trajectory(1:2,:));
angle = atan2(dir(3),dir(1));
th = [0 2*pi/3 4*pi/3]+angle;
r = 1;
pts0 = [chaotic_trajectory(1,1)+r*sin(th); chaotic_trajectory(1,3)+r*cos(th);];
dir = diff(chaotic_trajectory(end-1:end,:));
angle = atan2(dir(3),dir(1));
th = [0 2*pi/3 4*pi/3]+angle;
r = 1;
pts1 = [chaotic_trajectory(end,1)+r*sin(th); chaotic_trajectory(end,3)+r*cos(th);];

figure
setlatexlabels
imagesc(rho',XData=x,YData=z); 
set(gca,'ydir','normal');
colormap(flip(bone));
hold on
k=10;
plot(chaotic_trajectory(:,1),chaotic_trajectory(:,3),'r','LineWidth',lw);
fill(pts0(1,[1:end 1]),pts0(2,[1:end 1]),'r','EdgeColor','none');
fill(pts1(1,[1:end 1]),pts1(2,[1:end 1]),'r','EdgeColor','none');
xlabel('$x$','Interpreter','latex');
% ylabel('$z$','Interpreter','latex');
set(gcf,'color','w');
set(gca,'fontsize',fs);
axis equal;
ylim([0 50]);
xlim([-25 25]);
xticks([-20 0 20]);
yticks([]);

exportgraphics(gcf,'media/fig0b.pdf');

figure
setlatexlabels
imagesc(histogram',XData=x,YData=z); 
set(gca,'ydir','normal');
colormap(flip(bone));
hold on
% plot(chaotic_trajectory(:,1),chaotic_trajectory(:,3),'r','LineWidth',lw);
% fill(pts0(1,[1:end 1]),pts0(2,[1:end 1]),'r','EdgeColor','none');
% fill(pts1(1,[1:end 1]),pts1(2,[1:end 1]),'r','EdgeColor','none');
xlabel('$x$','Interpreter','latex');
% ylabel('$z$','Interpreter','latex');
set(gcf,'color','w');
set(gca,'fontsize',fs);
axis equal;
ylim([0 50]);
xlim([-25 25]);
xticks([-20 0 20]);
yticks([]);

exportgraphics(gcf,'media/fig0c.pdf');

end


