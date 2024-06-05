% author: Mauro Morini  
% last modified: 05.06.24
% This script should create different meshes with different parameters and
% circular resonators and have a section to calculate the time independent
% resonance frequencies
clc;clear;close all;

%% Global variables 
r = 0.1;
H = 0.02;
h = 0.005;
hGrad = 1.3;

%% No resonator, standard wave guide
% Initializations
% recRes = {[12, 0; 12.25, 5.5],[12, 6; 12.25, 10]};
% H = 0.05;
% h = 0.05;
% hGrad = 1.3;
% kappa1 = 0.001;
% 
% % Wave guide initialization
% W = WaveGuide;
% W = W.setMeshsize(H,h,hGrad);
% W = W.updateModel();
% W.plotMesh()
% W = W.assembleMatrices();
% W = W.lumpM();

%% 1 Circular resonator in square lattice
% Initializations
p0 = [0,0];
xNum = 3;
yNum = 3;
xLen = 1;
yLen = 1;
kappa1 = 0.001;

% Wave guide 
W = WaveGuide([p0; p0(1) + xNum*xLen, p0(2) + yNum*yLen]);
c0 = [p0(1)+xLen/2, p0(2)+yLen/2];
W = W.setMeshsize(H,h,hGrad);
for i = 1:xNum
    for j = 1:yNum
        W = W.addCircRes(c0+[(i-1)*xLen, (j-1)*yLen], r);
    end
end
W = W.updateModel();
W.plotMesh()
W = W.assembleMatrices();
W = W.lumpM();
calculateSafeResFreq(W,"square_lattice_grid_with_one_circ_res2.mat",r);

%% 2 Circular resonator in square lattice
% Initializations
p0 = [0,0];
xNum = 3;
yNum = 3;
xLen = 1;
yLen = 1;

% Wave guide 
W = WaveGuide([p0; p0(1) + xNum*xLen, p0(2) + yNum*yLen]);
c0 = [p0(1)+xLen/2 - 1.2*r, p0(2)+yLen/2];
c1 = [p0(1)+xLen/2 + 1.2*r, p0(2)+yLen/2];
W = W.setMeshsize(H,h,hGrad);
for i = 1:xNum
    for j = 1:yNum
        W = W.addCircRes(c0+[(i-1)*xLen, (j-1)*yLen], r);
        W = W.addCircRes(c1+[(i-1)*xLen, (j-1)*yLen], r);
    end
end
W = W.updateModel();
W.plotMesh()
W = W.assembleMatrices();
W = W.lumpM();
calculateSafeResFreq(W,"square lattice grid with two circ res2.mat",r);

%% Circular resonator in hexagonal lattice
% Initializations
p0 = [0,0];
xNum = 3;
yNum = xNum;
xLen = 3;
yLen = 2*sqrt(3);

% Wave guide 
W = WaveGuide([p0 - [0, yNum*yLen/2]; p0 + [xNum*xLen, yNum*yLen/2]]);

W = W.setMeshsize(H,h,hGrad);
for i = 1:xNum
    for j = 1:yNum
        c = p0 + (i-1)*[xLen, yLen]/2 + (j-1)*[xLen, -yLen]/2;
        W = W.addCircRes(c + [xLen/3, 0], r);
        W = W.addCircRes(c + [2*xLen/3, 0], r);
    end
end
W = W.updateModel();
W.plotMesh()
W = W.assembleMatrices();
W = W.lumpM();
calculateSafeResFreq(W,"2 resonator in hexagonal lattice grid2.mat",r);

%% Circular resonator in honeycomb lattice
% Initializations
p0 = [0,0];
xNum = 3;
yNum = xNum;
xLen = 3;
yLen = 2*sqrt(3);

% Wave guide 
W = WaveGuide([p0 - [0, yNum*yLen/2]; p0 + [xNum*xLen, yNum*yLen/2]]);

W = W.setMeshsize(H,h,hGrad);
for i = 1:xNum
    for j = 1:yNum
        c = p0 + (i-1)*[xLen, yLen]/2 + (j-1)*[xLen, -yLen]/2;
        c1 = c + [xLen/3 + 3*r, 0];
        c2 = c + [xLen/3,0] + 3*r*[cos(2*pi/3), sin(2*pi/3)];
        c3 = c + [xLen/3,0] + 3*r*[cos(4*pi/3), sin(4*pi/3)];
        c4 = c + 2*[xLen/3,0] + 3*r*[cos(pi/3), sin(pi/3)];
        c5 = c + 2*[xLen/3,0] - 3*r*[xLen/3, 0];
        c6 = c + 2*[xLen/3,0] + 3*r*[cos(5*pi/3), sin(5*pi/3)];
        W = W.addCircRes(c1, r);
        W = W.addCircRes(c2, r);
        W = W.addCircRes(c3, r);
        W = W.addCircRes(c4, r);
        W = W.addCircRes(c5, r);
        W = W.addCircRes(c6, r);
    end
end
W = W.updateModel();
W.plotMesh()
W = W.assembleMatrices();
W = W.lumpM();
calculateSafeResFreq(W,"honeycomb lattice grid2.mat",r);

%% Calculate time independent resonance frequencies
function calculateSafeResFreq(W,name,r)
    kappa1 = 0.001;
    [A,M] = W.getGlobMat([1, 1/kappa1], [1, 1/kappa1]);
    eigMax = sqrt(eigs(A,M, 10));
    eigMin = sqrt(eigs(A,M, 10, "smallestabs"));
    eigMid = sqrt(eigs(A,M, 10, eigMax(1)));
    save(name, "eigMax","eigMid","eigMin","W","r")
end

%% computation of all together
r = [0.1, 0.01];
meshUpdateParam = [0.2, 0.005, 1.3];
p0 = [0,0];
Num = [3,3];
Len = [1, 1; 3, 2*sqrt(3)];
idxR = [1,1,1,1,2,2,2,2];
idxLen = [1,1,2,2,1,1,2,2];
idxName = [1:4,1:4];
Name = ["1_in_Square", "2_in_Square", "2_in_Quadri", "honeycomb"];
kappa = 0.001;

for i = 1:8
    W = WaveGuide(p0, meshUpdateParam, Num, Len(idxLen(i),:),r(idxR(i)), Name(idxName(i)));
    W = W.assembleMatrices;
    W = W.lumpM;
    W.plotMesh;
    eigMax = W.computeResonanceFrequencies([1, 1/kappa], [1, 1/kappa], 10, "largestabs");
    eigMin = W.computeResonanceFrequencies([1, 1/kappa], [1, 1/kappa], 10, "smallestabs");
    eigMid = W.computeResonanceFrequencies([1, 1/kappa], [1, 1/kappa], 10, eigMax(1)/2);
    filename = Name(idxName(i)) + idxR(i) + ".mat";
    save(filename, "eigMax","eigMid","eigMin","W","r")
end