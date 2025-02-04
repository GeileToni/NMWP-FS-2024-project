% author: Mauro Morini  
% last modified: 09.05.24
function z = inShape2D(p, S)
% checks if given 2D point matrix is inside (or on the edge) of the given
% shape, which can either be a cluster of polygons or a cluster of non 
% intersecting circles
%
% Inputs:
% p: (n1,2) point matrix with x values in first and y in second column
% S: (n2,i) object matrix, for polygons i = 2, for circles i = 3 where the
%       first two columns are x and y coordinates of the nodes, for circles 
%       the last column denotes the radius, each object is separated by a 
%       [NaN,NaN] row for polygons, for cirlces each row is its unique
%       entity
%
% Outputs:
% z: logical (n1,1) index array of points in the shape

% Initializations
[n, k] = size(S);


% polygon case
if k == 2
    z = inpolygon(p(:,1), p(:,2),S(:,1), S(:,2));
    return;
end

% circle case
z = false(size(p,1),1);
for i = 1:n
    pNorm = vecnorm(p-S(i,1:2),2,2);
    z = z + (pNorm <= S(i,3));
end
z = logical(z);
end