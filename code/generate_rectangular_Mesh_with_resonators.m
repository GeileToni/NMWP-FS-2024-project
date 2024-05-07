% author: Mauro Morini
% last modified 14.04.24
clc;clear;close all;

R1coord = [12, 0; 12.25, 0; 12.25, 5.5; 12, 5.5];
R2coord = [12, 6; 12.25, 6; 12.25, 10; 12, 10];
RTotCoord = [R1coord; NaN, NaN; R2coord];
pgon = polyshape(RTotCoord(:,1),RTotCoord(:,2));
dl = gen_mesh(1);
model = createpde;
geometryFromEdges(model, dl);
generateMesh(model, GeometricOrder='linear', ...
    Hface={[2,3],0.03}, Hmax=0.5, Hgrad=1.3);

% get coordinate, edge and connectivity matrix
[p,e,t] = meshToPet(model.Mesh);
% p = p';
% t = t';
% t = t(:, 1:3);
% e = e';
% e = e(:, 1:2);

figure(2)
pdeplot(model)
hold on
% h = plot(pgon);
% h.FaceColor = 'green';
fill(R1coord(:,1),R1coord(:,2),'g')
fill(R2coord(:,1),R2coord(:,2),'g')
axis equal;axis tight;
hold off 