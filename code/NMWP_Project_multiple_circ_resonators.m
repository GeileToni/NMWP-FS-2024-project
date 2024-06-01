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
Ratio = 1.5;
Lmin = min(Lx,Ly);
xLengthRes = Lx/(2*xNum + 1);
yLengthRes = Ly/(2*yNum + 1);
xLengthSpace = xLengthRes*(1/Ratio);
yLengthSpace = yLengthRes*(1/Ratio);
rRes = min(xLengthRes,yLengthRes)*Ratio/2;
xLengthRes = xLengthRes*Ratio;
yLengthRes = yLengthRes*Ratio;


backgNodes = [0, 0; Lx, 0; Lx, Ly; 0, Ly];
resNodeLoc = [rRes, rRes, rRes]';

% create resNode for mesh creation
resNode = zeros(3, xNum*yNum);
resNode(3,:) = rRes;
counter = 1;
for i = 1:yNum
    for j = 1:xNum
        resNode(1,counter) = resNodeLoc(1) + j*xLengthSpace + (j-1)*xLengthRes;
        resNode(2,counter) = resNodeLoc(2) + i*yLengthSpace + (i-1)*yLengthRes;
        counter = counter + 1;        
    end
end

% possibly take out a resonator
resNode = [resNode(:,1:5),resNode(:,7:end)];

% create mesh
[p,t,e, model] = createCircleResonatorMesh(backgNodes, resNode,Hmax,hres);

pdeplot(model)

% create ResTot for inShape
RTotCoord = resNode';

% functions
kappa = @(x,y) 1*~inShape2D([x, y], RTotCoord) + ...
        0.2*inShape2D([x, y], RTotCoord);

roh_r = @(t) 0.1/(1+0.2*cos(pi/2*t));
%roh_r = @(t) 0.2*ones(size(t));
roh = @(x,y,t) 1*~inShape2D([x, y], RTotCoord) + ...
        roh_r(t)*inShape2D([x, y], RTotCoord);
g = @(x,y,t) sin(omega*(x - t));
u0 = @(x,y) 0*x;
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
