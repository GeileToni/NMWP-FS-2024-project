% author: Mauro Morini  
% last modified: 09.05.24
function [uh_T, tvec] = waveEqLF2D(u0, u1, g, T, dt, p, e, t, kappa, roh, RTotCoord)
% solves the wave equation in 2D with Leap-Frog time discretization and
% mass lumping and FEM spacial discretization and returns u(x, T) at end
% time T with 4 boundary sections on which  2 have homogeneous neumann
% b.c. one has non-homogeneous neumann condition and one has an absorbing
% b.c. of first order (du/dt = -cdu/dn) on a rectangular mesh
% PDE : (1/kappa d^2/dt^2 - grad_x 1/roh grad_x)u(x,t) = 0
%
% Inputs :
% u0 : function handle for dirichlet initial condition u(x,0) = u0
% u1 : function handle for neumann initial condition du/dt (x,0) = u1
% g : function handle for neumann-b.c.
% T : scalar final time T, can be given as a row vector with times which should
%       be recorded (has to be sorted and actual recorded time t_actual 
%       will be different from requested time by an error of dt)
% dt : scalar time step size
% p : nPx2 coordinate matrix with points in rows
% t : nEx3 connectivity matrix with elements in rows
% e : (nE,2) edge matrix containing two connecting edge point indices in
%       each row, of the boundary edges of the domain
% kappa : function handle kappa(x,y) for time dependant part
% roh : function handle roh(x,y,t) for grad part
%       RTotCoord: (nRp, 2) matrix with coordinates of the vertices 
%       which contain the resonators, each resonator 
%       is seperated by a [NaN,NaN] row (assumes resonators to be polygonal)
%
% Outputs :
% uh_T : (nP,1) numerical solution at time T, if T is a vector then uh_T
%           has size (nP, length(T)) 

% Initializations
hasResonators = false;
RA_timedep = true;
TisVec = false;
Lx = max(p(:, 1));                % x length of wave guide

% case if T is a vector
if size(T, 2) > 1
    TisVec = true;
    TOld = T;
    iterT = 1;
    T = T(end);
    uh_T = zeros(size(p, 1), length(TOld));

    % if TOld(1) < 2*dt
    %     error("smallest time in vector T is too small: Tmin = " + Tsorted(1))
    % end
end

if ~exist('kappa','var') 
    kappa = @(x,y) ones(size(x)); 
end
if ~exist('roh','var') 
    roh = @(x,y,t) ones(size(x)); 
    RA_timedep = false;
end
c = @(x,y,t) sqrt(kappa(x,y)./roh(x,y,t));

if exist('RTotCoord','var')
    hasResonators = true;
end

% find boundary sections        
pe1x = p(e(:, 1), 1);
pe1y = p(e(:, 1), 2);
pe2x = p(e(:, 2), 1);
pe2y = p(e(:, 2), 2);

% gamma1 x = 0
gamma1Idx = logical((pe1x == 0).*(pe2x == 0));
gamma1 = e(gamma1Idx, :);

% gamma4 x = Lx
gamma2Idx = logical((pe1x == Lx).*(pe2x == Lx));
gamma2 = e(gamma2Idx, :);

% find elements t(i,:) which are in the resonator 
if hasResonators
    pt1 = p(t(:,1),:);
    pt2 = p(t(:,2),:);
    pt3 = p(t(:,3),:);
    % resIdx = logical( ...
    %     inpolygon(pt1(:,1),pt1(:,2),RTotCoord(:, 1), RTotCoord(:, 2)).*...
    %     inpolygon(pt2(:,1),pt2(:,2),RTotCoord(:, 1), RTotCoord(:, 2)).*...
    %     inpolygon(pt3(:,1),pt3(:,2),RTotCoord(:, 1), RTotCoord(:, 2)));
    resIdx = logical( ...
        inShape2D(pt1,RTotCoord).*...
        inShape2D(pt2,RTotCoord).*inShape2D(pt3,RTotCoord));
    xRes = RTotCoord(1,1);
    yRes = RTotCoord(1,2);
else
    resIdx = false(size(t,1),1);
    xRes = 0;
    yRes = 0;
end
% Assemble stiffness matrix inside resonator and outside
Ares = FEM2D.stiffnessMatrix2D(p, t(resIdx,:));
Abackg = FEM2D.stiffnessMatrix2D(p, t(~resIdx,:));
A = Ares/roh(xRes,yRes,0) + Abackg;

% % test
% Astand = stiffnessMatrix2D(p, t, @(x,y) 1/roh(x,y,0));
% norm(full(Astand - A),'fro')

% Assemble and lump remaining matrices
M = FEM2D.massMatrix2D(p, t, @(x,y) 1/kappa(x,y));
R = FEM2D.neumannMassMatrix2D(p, gamma2,@(x,y) ((c(x,y,0).*roh(x,y,0)).^(-1)));
G = FEM2D.neumannLoadVector2D(p, gamma1, @(x,y) g(x,y,0)./roh(x,y,0));
MLump = diag(sum(M, 2));
RLump = diag(sum(R, 2));
M = MLump;
R = RLump;

% calculate stable dt
lMax = eigs(A,1);
lMin = eigs(M,1,0);
dt = sqrt(lMin/lMax);
tvec = 0:dt:T;
N = length(tvec) - 1;

% initial conditions
uhPrev = u0(p(:, 1), p(:, 2));
uhNow = (2/dt^2*M)\(G - (A - 2/dt^2*M)*uhPrev + (1/dt^2*M - 1/(2*dt)*R)*2*dt*u1(p(:, 1), p(:, 2)));

for i = 1:N

    % construct system
    if RA_timedep
        % not used since roh is just 1 on the boundary
        % R = neumannMassMatrix2D(p, gamma2,@(x,y)((c(x,y,dt*i).*roh(x,y,dt*i)).^(-1)));   
        % A = stiffnessMatrix2D(p, t, @(x,y) 1/roh(x,y,dt*i));
        A = Ares/roh(xRes,yRes,dt*i) + Abackg;
    end
    G = FEM2D.neumannLoadVector2D(p, gamma1, @(x,y) g(x,y,dt*i)./roh(x,y,dt*i));
    LHS = (1/dt^2*M + 1/(2*dt)*R);
    RHS = G - (A - 2/dt^2*M)*uhNow - (1/dt^2*M - 1/(2*dt)*R)*uhPrev;

    % solve system
    uhNext = LHS\RHS;
    uhPrev = uhNow;

    if TisVec && abs(dt*i - TOld(iterT)) < (dt - eps)
        uh_T(:, iterT) = uhNow;
        iterT = iterT + 1;
        disp("Computation for t = " + dt*i + " completed")
    end
    uhNow = uhNext;
end

if ~TisVec
    uh_T = uhNext;
end

if nargout > 1
    tvec = 0:dt:T;
end

% last step with smaller stepsize to arrive at T !!!!!NOT GOOD!!!!!
% if N*dt < T 
%     dt2 = T - N*dt;
% 
%     % construct system
%     if RA_timedep
%         % R = neumannMassMatrix2D(p, gamma2,@(x,y) ((c(x,y,T).*roh(x,y,T)).^(-1)));
%         A = stiffnessMatrix2D(p, t, @(x,y) 1/roh(x,y,T));
%     end
%     G = neumannLoadVector2D(p, gamma1, @(x,y) g(x,y,T)./roh(x,y,T));
%     LHS = (1/dt2^2*M + 1/(2*dt2)*R);
%     RHS = G - (A - 2/dt2^2*M)*uhNow - (1/dt2^2*M - 1/(2*dt2)*R)*uhPrev;
% 
%     % solve system
%     uhNext = LHS\RHS;
%     if ~TisVec
%         uh_T = uhNext;
%     else
%         uh_T(:,iterT) = uhNext;
%         if iterT ~= length(TOld)
%             error("iterT is not the end time, something went wrong: iter: " + iter + " .. T: " +T)
%         end
%     end
% 
% end
end