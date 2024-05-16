% author: Mauro Morini  
% last modified: 08.05.24
function [p,t,e, model] = createCircleResonatorMesh(backgNodes, resNodes,Hmax,hres)
% creates a mesh of polygonal shape with a background medium and multiple
% resonators inside, all resonators are circles 
%
% Inputs:
% backgNodes: (n1,2) containing nodes of background medium 
%            (i.e. the outlines of the mesh) in first column the x in the
%            second column the y values
% resNodes: (3, n2) containing nodes of the resonators as circles,
%               each column is first the x then the y values of the center
%               and finally the radius (which is the same)
% Hmax: scalar maximal mesh size for the background medium
% hres: scalar meshsize inside of the resonators
%
% Outputs:
% standard p,e,t and model

backg = [3;4;backgNodes(:,1); backgNodes(:,2)];
res = [ones(1,size(resNodes,2)); resNodes;zeros(6,size(resNodes,2))];
gd = [backg, res];
dl = decsg(gd);

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