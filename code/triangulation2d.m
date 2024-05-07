% author: Mauro Morini
% last modified: 06.11.23
function [p, t, e] = triangulation2d(p)
% creates a 2d triangulation with given points p
% Outputs : 
% p :   coordinate matrix nPx2 containing points (x,y) in rows, 
%       representing the verteces of the triangles  
% t :   connectivity matrix nTx3 each row representing one element 
%       and each column representing the local numbering (0,1,2) and the 
%       entry the global numbering
% e :   connectivity matrix nEx2 for edges on the boundary
%
% Input : 
% p :   same as output but with possible duplicates
DT = delaunayTriangulation(p);

% eliminate duplicates
p = DT.Points;
t = DT.ConnectivityList;
e = freeBoundary(DT);

end