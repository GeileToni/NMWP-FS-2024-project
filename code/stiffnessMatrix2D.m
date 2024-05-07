% author: Mauro Morini
% last modified 24.04.24
function A = stiffnessMatrix2D(p, t, c)
% calculates global stiffness matrix for linear FEM in 2D int_omega
% a(x)*phi_i*phi_j dx
% Inputs : 
% p : nPx2 coordinate matrix with points in rows
% t : nEx3 connectivity matrix with elements in rows
% c : if exists function handle of bilinear form
% Outputs : 
% A : global stiffness matrix nPxnP

if ~exist('c','var') 
    c = @(x,y) 1; 
end

nP = size(p, 1);
nE = size(t, 1);
A = sparse(nP, nP);

% iterate over elements
for k = 1:nE
    % Element and points
    K = t(k, :);
    p0 = p(K(1), :).';
    p1 = p(K(2), :).';
    p2 = p(K(3), :).';

    cent = mean([p0, p1, p2], 2);       % center point

    % calculate J_K elementwise jacobian of the transformation
    Jk = [p1-p0, p2-p0];

    K_area = abs(det(Jk))/2;
    DN = [-1, -1; 1, 0; 0, 1];

    % elementwise local stiffness matrix
    A_loc = c(cent(1),cent(2))*K_area*DN*inv(Jk'*Jk)*DN';

    % Assembling
    A(t(k, :), t(k, :)) = A(t(k, :), t(k, :)) + A_loc;
end
end