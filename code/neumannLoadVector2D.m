% author: Mauro Morini (adapted from "The Finite Element Theory" by
% Larson/Bengzon chapter 4.6)
% last modified 16.04.24
function L = neumannLoadVector2D(p, e, f)
% assembles load vector for neumann b.c. : du/dn = f on boundary Gamma given by
% the edge matrix e for 2d linear FE. In the weak formulation, the above
% b.c. has the following impact on the LHS: int_Gamma v*f ds for a
% testfunction v. 
%
% Inputs:
% p : (nP,2) point matrix with x values in the first and y values in the
%       second column
% e : (nE,2) edge matrix containing two connecting edge point indices in
%       each row, of the boundary edges of the domain
% f : function handle of neumann b.c.
%
% Outputs:
% L : (nP,1) neumann load vector containing the pointwise impact of the
%       b.c. in each entry (zero for internal points)

% Initializations
nP = size(p, 1);
nE = size(e, 1);
L = zeros(nP,1);

for i = 1:nE
    E = e(i,:);
    x1 = p(E(1),:);                               % edge point 1
    x2 = p(E(2),:);                               % edge point 2
    e_h = norm(x1 - x2, 2);                       % edge length
    cent = mean([x1', x2'], 2);

    % local element contributions
    Lloc = e_h/2*f(cent(1),cent(2))*[1;1];          
    L(E) = L(E) + Lloc;
end
end