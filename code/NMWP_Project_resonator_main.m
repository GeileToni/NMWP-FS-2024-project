% author: Mauro Morini  
% last modified: 01.06.24
clc;clear;close all;

% Initializations
recRes = {[12, 0; 12.25, 5.5],[12, 6; 12.25, 10]};
T = 0.5:0.5:30;
omega = 2*pi;
H = 0.1;
h = 0.05;
hGrad = 1.3;

% functions
rho_r = @(t) 0.1/(1+0.2*cos(pi/2*t));
kappa1 = 0.001;
g = @(x,y,t) sin(omega*(x - t));
u0 = @(x,y) x*0;
u1 = u0;

% Wave guide initialization
W = WaveGuide;
W = W.setMeshsize(H,h,hGrad);
% for i = 1:length(recRes)
%     W = W.addRectRes(recRes{i});
% end
W = W.addCircRes([W.bbox(2,1), W.bbox(2,2)]/2,0.1);
W = W.updateModel();
W.plotMesh()
W = W.assembleMatrices();
W = W.lumpM();
[p,e,t] = W.getPet;

% calculate natural frequencies 
[A,M] = W.getGlobMat([1, 1/kappa1], [1, 1/kappa1]);
dTot = -sqrt(eigs(A,M));
resIdx = W.getPointsInRes;
dRes = -sqrt(eigs(W.Ares(resIdx,resIdx)/kappa1,W.Mres(resIdx,resIdx)/kappa1));
%% calculate time dependant solution uh
g = @(x,y,t) sin(dRes(1)*(x - t));
% find boundary sections        
pe1x = p(e(:, 1), 1);
pe1y = p(e(:, 1), 2);
pe2x = p(e(:, 2), 1);
pe2y = p(e(:, 2), 2);

% gamma1 x = 0
gamma1Idx = logical((pe1x == 0).*(pe2x == 0));
gamma1 = e(gamma1Idx, :);

% gamma4 x = Lx
gamma2Idx = logical((pe1x == W.bbox(2,1)).*(pe2x == W.bbox(2,1)));
gamma2 = e(gamma2Idx, :);

% matrices
[A,M] = W.getGlobMat([1, 1/rho_r(0)], [1, 1/kappa1]);
R = FEM2D.neumannMassMatrix2D(p, gamma2);
R = diag(sum(R, 2));            % lumped R
G = FEM2D.neumannLoadVector2D(p, gamma1, @(x,y) g(x,y,0));

% calculate stable dt
lMax = eigs(A,1);
lMin = eigs(M,1,0);
dt = 0.9*2*sqrt(lMin/lMax);

% initial conditions
uhPrev = u0(p(:, 1), p(:, 2));
uhNow = (2/dt^2*M)\(G - (A - 2/dt^2*M)*uhPrev + (1/dt^2*M - 1/(2*dt)*R)*2*dt*u1(p(:, 1), p(:, 2)));

% Leap Frog
i = 1;
iterT = 1;
uh_T = zeros(size(uhNow));
while (i+1)*dt <= T(end)
    A = W.getGlobStiffness([1, 1/rho_r(i*dt)]);
    G = FEM2D.neumannLoadVector2D(p, gamma1, @(x,y) g(x,y,dt*i));
    LHS = (1/dt^2*M + 1/(2*dt)*R);
    RHS = G - (A - 2/dt^2*M)*uhNow - (1/dt^2*M - 1/(2*dt)*R)*uhPrev;

    % solve system
    uhNext = LHS\RHS;

    % update
    uhPrev = uhNow;
    i = i+1;
    uhNow = uhNext;
    
    % save uhNow if close enough to wanted T value
    if abs(dt*i - T(iterT)) < (dt - eps)
        uh_T(:,iterT) = uhNow;
        T(iterT) = dt*i;
        iterT = iterT+1;
    end
end