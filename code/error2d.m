% author: Mauro Morini  
% last modified: 18.11.23
function L2err = error2d(p, t, u, uh)
% calculates the L^2 error of a linear FE soltion in 2D
%
% Inputs : 
% p : nPx2 coordinate matrix with points in rows
% t : nEx3 connectivity matrix with elements in rows
% u : function handle @(x,y) of exact solution u
% uh : nPx1 column vector of FE solution
%
% Output : 
% L2err : scalar value of L^2 error (||u-uh||_{L^2(omega)})

% Initalizations
[nE, r] = size(t);
nP = size(p, 1);
L2err = 0;
KhatArea = 1/2;

% weights of QF
w = (1/3)*[1, 1, 1];

% nodes of QF 2xr
xi = [0, 0; 1, 0; 0, 1]';

% shape functions
N = {@(xi) 1 - [1,1]*xi, @(xi) [1,0]*xi, @(xi) [0,1]*xi};

% iterate over elements
for i = 1:nE
    
    % Element and points
    K = t(i, :);
    p0 = p(K(1), :).';
    p1 = p(K(2), :).';
    p2 = p(K(3), :).';

    % calculate J_K elementwise jacobian of the transformation
    Jk = [p1-p0, p2-p0];

    % element map
    Fk = @(xi) p0 + Jk*xi;

    KArea = abs(det(Jk))/2;

    % calculate Quadrature
    Qval = zeros(1, r);
    for j = 1:r
        uhVal = 0;
        for q = 1:r
            uhVal = uhVal + uh(t(i, q))*N{q}(xi(:, j));
        end
        x = Fk(xi(:, j));
        Qval(j) = abs(u(x(1), x(2)) - uhVal)^2;
    end
    L2err = L2err + (KhatArea*KArea)*Qval*w';
end
L2err = sqrt(L2err);
end