% author: Mauro Morini  
% last modified: 23.04.24
clc;clear;close all;
%%
% Initializations
h = 2^(-4)/2;
dt = 0.95*h^2;
T = 0.5;

% functions
u0 = @(x,y) x*0;
u1 = @(x,y) sin(pi*x);
g = @(x,y,t) sin(pi*(x - t));
u = @(x,y,t) sin(pi*x).*sin(pi*t)/pi;

% mesh
[p, t, e] = generateMesh2dUnitSquare(h);

uh_T = waveEqLF2D(u0, u1, g, T, dt, p, e, t);


figure(1);
tld = tiledlayout('flow');
nexttile
trisurf(t,p(:,1),p(:,2),uh_T, 'EdgeColor', 'none')
title("uh_T")

nexttile
trisurf(t,p(:,1),p(:,2),u(p(:,1),p(:,2), T), 'EdgeColor', 'none')
title("uh")

%% inpolygon tests
R1coord = [12, 0; 12.25, 0; 12.25, 5.5; 12, 5.5];
R2coord = [12, 6; 12.25, 6; 12.25, 10; 12, 10];
RTotCoord = [R1coord; NaN, NaN; R2coord];
% 
% kappa = @(x,y) inpolygon(x,y,R1coord(:, 1), R1coord(:, 2)) ...
%         + inpolygon(x,y,R2coord(:, 1), R2coord(:, 2));

kappa = @(x,y) inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));

xq = rand(25000, 1)*20; yq = rand(25000, 1)*10;
idx = ~kappa(xq,yq);
figure(2);
plot(xq(idx),yq(idx), 'x', [R1coord(:,1);R2coord(:,1)], [R1coord(:,2);R2coord(:,2)], 'o')

kappa = @(x,y) 1*~inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2)) + ...
        0.2*inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));
roh_r = @(t) 0.1/(1+0.2*cos(pi/2*t));
roh = @(x,y,t) 1*~inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2)) + ...
        roh_r(t)*inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));

[X,Y] = meshgrid(linspace(0, 20, 2000), linspace(0, 10, 1000));
V = roh(X,Y,0.5);
figure(3);
surf(X,Y,V, 'LineStyle','none')
view(2)

figure(4)
pgon = polyshape(RTotCoord(:,1),RTotCoord(:,2));
plot(pgon)
% xlim([0,20])
% ylim([0,10])