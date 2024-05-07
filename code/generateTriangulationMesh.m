% author: Mauro Morini
% last modified 05.12.23
function [p,t,e, model] = generateTriangulationMesh(hmax, gd, sf, ns, refine)
% creates a mesh and the p,e,t matrices from given geometry inputs
dl = decsg(gd, sf, ns);

% create a model object and add geom to it then generate mesh for a single
% maximal meshsize h
model = createpde;
geometryFromEdges(model, dl);
generateMesh(model, "GeometricOrder","linear", "Hmax",hmax);

% refine mesh
if exist('refine', 'var')  
    model = refinePDEMmesh(model);
end
% get coordinate, edge and connectivity matrix
[p,e,t] = meshToPet(model.Mesh);
p = p';
t = t';
t = t(:, 1:3);
e = e';
e = e(:, 1:2);
end