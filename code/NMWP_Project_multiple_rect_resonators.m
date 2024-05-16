% author: Mauro Morini  
% last modified: 08.05.24
clc;clear;close all;

% Initialization 
Hmax = 0.1;
hres = 0.1;
r = 0.4;
dt = 0.005;
omega = 2*pi;
T = 0.5:0.2:30;

% mesh nodes
Lx = 20;
Ly = 10;
xNum = 4;
yNum = 3;
xRatio = 1.5;
yRatio = 1.5;
xLengthRes = Lx/(2*xNum + 1);
yLengthRes = Ly/(2*yNum + 1);
xLengthSpace = xLengthRes*(1/xRatio);
yLengthSpace = yLengthRes*(1/yRatio);
xLengthRes = xLengthRes*xRatio;
yLengthRes = yLengthRes*yRatio;

backgNodes = [0, 0; Lx, 0; Lx, Ly; 0, Ly];
resNodeLoc = [0, 0; xLengthRes, 0; xLengthRes, yLengthRes; 0, yLengthRes];

% create resNode for mesh creation
resNode = zeros(8, xNum*yNum);
counter = 1;
for i = 1:yNum
    for j = 1:xNum
        resNode(1:4,counter) = resNodeLoc(:,1) + j*xLengthSpace + (j-1)*xLengthRes;
        resNode(5:8,counter) = resNodeLoc(:,2) + i*yLengthSpace + (i-1)*yLengthRes;
        counter = counter + 1;        
    end
end

% possibly take out a resonator
resNode = [resNode(:,1:5),resNode(:,7:end)];

% create mesh
[p,t,e, model] = createRectangleResonatorMesh(backgNodes, resNode,Hmax,hres);

%pdeplot(model)

% create ResTot for inpolygon
counter = 1;
nR = size(resNode,2);
for i = 1:nR
    RTotCoord(counter:(counter+3),1) = resNode(1:4,i);
    RTotCoord(counter:(counter+3),2) = resNode(5:8,i);
    RTotCoord(counter+4,:) = [NaN, NaN];
    counter = counter+5;
end
% pgon = polyshape(ResTot);
% plot(pgon)

% functions
kappa = @(x,y) 1*~inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2)) + ...
        0.2*inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));

roh_r = @(t) 0.1/(1+0.2*cos(pi/2*t));
roh_r = @(t) 0.2*ones(size(t));
roh = @(x,y,t) 1*~inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2)) + ...
        roh_r(t)*inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));
g = @(x,y,t) sin(omega*(x - t));
u0 = @(x,y) x*0;
u1 = u0;

%% compute numerical solution
uh_T = waveEqLF2D(u0, u1, g, T, dt, p, e, t, kappa, roh,RTotCoord);

%% plot
for i = length(T):-1:1
    figure;
    tld = tiledlayout('flow');
    nexttile
    trisurf(t,p(:,1),p(:,2),uh_T(:,i), 'EdgeColor', 'none')
    title("uh_T, t = " + T(i))
    %view(2)
end
