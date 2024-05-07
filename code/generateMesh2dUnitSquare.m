% author: Mauro Morini
% last modified: 06.11.23
function [p, t, e] = generateMesh2dUnitSquare(h)
% creates an equidistant mesh and a triangularization on the unit square
% (0,1)x(0,1)
%
% Outputs : 
% p :   coordinate matrix nPx2 containing points (x,y) in rows, 
%       representing the verteces of the triangles  
% t :   connectivity matrix nTx3 each row representing one element 
%       and each column representing the local numbering (0,1,2) and the 
%       entry the global numbering
% e :   connectivity matrix nEx2 for edges on the boundary of the form
%       (i,j), where i and j are indexes of points on the boundary
%
% Inputs : 
% h :   meshsize (equidistant)

[X, Y] = meshgrid(0:h:1, 1:-h:0);
N = size(X, 1);

% construct points p
p = zeros(N^2, 2);
k = 1;
for i = 1:N
    for j = 1:N
        % x coordinate
        p(k, 1) = X(i,j);

        % y coordinate
        p(k, 2) = Y(i,j);
        k = k + 1;
    end
end

[p, t, e] = triangulation2d(p);
end