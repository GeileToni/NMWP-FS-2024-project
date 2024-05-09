% author: Mauro Morini  
% last modified: 08.05.24
function [p,t,e, model] = createRectangleResonatorMesh(backgNodes, resNodes,Hmax,hres)
% creates a mesh of polygonal shape with a background medium and multiple
% resonators inside, all resonators are rectangular and are contained
% in resNodes, they are separated by NaN
%
% Inputs:
% backgNodes: (n1,2) containing nodes of background medium 
%            (i.e. the outlines of the mesh) in first column the x in the
%            second column the y values
% resNodes: (8, n2) containing nodes of the resonators as rectangles, each 
%               each column is first the 4 x then the 4 y values
% Hmax: scalar maximal mesh size for the background medium
% hres: scalar meshsize inside of the resonators
%
% Outputs:
% 

backg = [backgNodes(:,1); backgNodes(:,2)];
gd = [backg, resNodes];
gd = [3*ones(1,size(gd,2));4*ones(1,size(gd,2));gd];
% sf = 'all';
% ns = num2cell(1:size(gd,2));
% %ns = num2str((1:size(gd,2)).', '%d')';
% ns = string(1:size(gd,2));
dl = decsg(gd);
% dl = decsg(gd,sf,ns);

figure(1)
pdegplot(dl,"FaceLabels","on")

model = createpde;
geometryFromEdges(model, dl);
% create vector without face number of background medium
bgFace = 4;
resFaceIdx = [1:(bgFace-1), (bgFace+1):size(gd,2)];
generateMesh(model, "GeometricOrder","linear", Hmax=Hmax,Hface={resFaceIdx, hres},Hgrad=1.3);
% generateMesh(model, "GeometricOrder","linear", Hmax=Hmax);

[p,e,t] = meshToPet(model.Mesh);
p = p';
t = t';
t = t(:, 1:3);
e = e';
e = e(:, 1:2);

end