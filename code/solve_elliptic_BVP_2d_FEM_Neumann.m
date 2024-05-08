% author: Mauro Morini
% last modified 14.04.24
function uh = solve_elliptic_BVP_2d_FEM_Neumann(p, t, e, f, g)
% solves the problem -laplace(u) + u = f on the unit square in 2D with
% neumann b.c du/dn = g using linear FEM
%
% Inputs : 
% p : nPx2 coordinate matrix with points in rows
% t : nEx3 connectivity matrix with elements in rows
% f : function handle RHS
% g : function handle, g:omega->R, for boundary condition
% e : connectivity matrix for the free edges (edges on the boundary)
%
% Outputs : 
% uh : nPx1 coefficient vector of the FEM solution

% Initializations
nP = size(p, 1);

% assemble matrices
A = FEM2D.stiffnessMatrix2D(p,t);
M = FEM2D.massMatrix2D(p,t);
L = FEM2D.loadVector2D(p,t,f);
G = FEM2D.neumannLoadVector2D(p,e,g);

% % find interior points of triangularization 
% idxExt = reshape(e, [], 1);
% idxExt = ismember(1:nP, idxExt);
% idxInt =~ idxExt;
% 
% % incorporate b.c.
% G(idxInt) = 0;

% solve system
uh = (A + M)\(L + G);
end