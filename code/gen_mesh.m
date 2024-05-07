% author: Mauro Morini
% last modified 14.04.24
function geom = gen_mesh(flag)
% 
%   Inputs :
%   flag : either 0 or 1 deciding if figure should be plotted or not
% 
%   Outputs : 
%   geom : decomposed geometriy matrix for reference see decsg
%   documentation

% meshdefining points x-values, y-values respectively
x = [-1,1,1,0,0,-1];
y = [-1,-1,1,1,0,0];
% components for decsg (can easily be acquired using the modeler app/
% pdepoly(x,y))
% gd = [2,6,-1,1,1,0,0,-1,-1,-1,1,1,0,0]';
% ns = [80;49];
% sf = 'P1';
x1 = [0, 20, 20, 0];
y1 = [0, 0, 10, 10];
r1x = [12, 12.25, 12.25, 12];
r1y = [0, 0, 5.5, 5.5];
r2x = [12, 12.25, 12.25, 12];
r2y = [6, 6, 10, 10];
gd =   [3, 4, x1, y1;
        3, 4, r1x, r1y;
        3, 4, r2x, r2y]';
ns = ["F1", "F2", "F3"];
sf = 'F1+F2+F3';
geom = decsg(gd,sf,ns);

if flag ~= 0 && flag ~= 1
    error('flag parameter should be either 0 or 1 \n%s', "flag: " + flag)
elseif flag == 0
    return
end

figure()
pdegplot(geom, 'EdgeLabels','on',"FaceLabels","on")
xlabel("x")
ylabel("y")
title("Resonator mesh")
end