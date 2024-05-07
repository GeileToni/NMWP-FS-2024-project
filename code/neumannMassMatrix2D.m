% author: Mauro Morini (adapted from "The Finite Element Theory" by
% Larson/Bengzon chapter 4.6)
% last modified 23.04.24
function R = neumannMassMatrix2D(p, e, c)
% assembles mass matrix for neumann b.c. on the edges int_gamma
% c*phi_i*phi_j ds
%
% Inputs:
% p : (nP,2) point matrix with x values in the first and y values in the
%       second column
% e : (nE,2) edge matrix containing two connecting edge point indices in
%       each row, of the boundary edges of the domain
% c : function handle from weak formulation default is 1
%
% Outputs:
% R :  (nP,nP) Neumann mass matrix with R_ij = int_{e_ij} c phi_i phi_j ds

if ~exist('c','var') 
    c = @(x,y) 1; 
end

% Initializations
nP = size(p, 1);
nE = size(e, 1);
R = sparse(nP,nP);

for i = 1:nE
    E = e(i,:);
    x1 = p(E(1),:);                               % edge point 1
    x2 = p(E(2),:);                               % edge point 2
    e_h = norm(x1 - x2, 2);                       % edge length
    cent = mean([x1', x2'], 2);

    % local element contributions, same as for mass matrix but here only
    % 2x2 because the integral is over e
    Rloc = e_h/6*[2, 1; 1 2]*c(cent(1), cent(2));          
    R(E, E) = R(E, E) + Rloc;
end
end