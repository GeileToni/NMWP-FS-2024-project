% author: Mauro Morini  
% last modified: 03.05.24
clc;clear;close all;

% Initialization 
H = 0.1;
h = 0.03;
r = 0.4;
dt = h^2*r;
dt = 0.005;
omega = 2*pi;
T = 0.5:0.2:30;

% create mesh
R1coord = [12, 0; 12.25, 0; 12.25, 5.5; 12, 5.5];
R2coord = [12, 6; 12.25, 6; 12.25, 10; 12, 10];
RTotCoord = [R1coord; NaN, NaN; R2coord];
dl = gen_mesh(0);
model = createpde;
geometryFromEdges(model, dl);
generateMesh(model, GeometricOrder='linear', Hface={[2,3],h}, ...
    Hmax=H, Hgrad=1.3);
[p,e,t] = meshToPet(model.Mesh);
p = p';
t = t';
t = t(:, 1:3);
e = e';
e = e(:, 1:2);

% functions
kappa = @(x,y) 1*~inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2)) + ...
        0.2*inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));
roh_r = @(t) 0.1/(1+0.2*cos(pi/2*t));
roh = @(x,y,t) 1*~inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2)) + ...
        roh_r(t)*inpolygon(x,y,RTotCoord(:, 1), RTotCoord(:, 2));
g = @(x,y,t) sin(omega*(x - t));
u0 = @(x,y) x*0;
u1 = u0;

%%
% uh_T = zeros(size(p, 1), length(T));
% for i = 1:length(T)
%     uh_T(:, i) = waveEqLF2D(u0, u1, g, T(i), dt, p, e, t, kappa, roh);
%     disp("Computation for t = " + T(i) + " completed")
% end
uh_T = waveEqLF2D(u0, u1, g, T, dt, p, e, t, kappa, roh,RTotCoord);

%% plot
for i = length(T):-1:1
    figure;
    tld = tiledlayout('flow');
    nexttile
    trisurf(t,p(:,1),p(:,2),uh_T(:,i), 'EdgeColor', 'none')
    title("uh_T, t = " + T(i))
    %view(2)
end
% figure(1);
% tld = tiledlayout('flow');
% Tloc = [59, 63, 67, 70, 75, 110];
% Tloc = 90;
% for i = Tloc
%     nexttile
%     figure
%     trisurf(t,p(:,1),p(:,2),uh_T(:,i), 'EdgeColor', 'none')
%     xlabel("x")
%     ylabel("y")
%     zlabel("uh")
%     title("uh_T, t = " + T(i))
%     colorbar;
%     axis equal; axis off; axis tight ; %colormap ('jet');
%     %view(2)
%     %view(270,90)
% end


% save("NMWP_Project_E4_workspace_backup_v1.06.mat")
% save("uh_T_E4_v1.06.mat","uh_T")
% save("NMWP_Project_resonator_mesh_data_v1.06.mat")

%% new adaptation using WaveGuide
clc;clear;close all;

% Initializations
R1coord = [12, 0; 12.25, 5.5];
R2coord = [12, 6; 12.25, 10];

% Wave guide initialization
W = WaveGuide();
W = W.addRectRes(R1coord);
W = W.addRectRes(R2coord);
W = W.updateModel();
W.plotMesh()
W = W.assembleMatrices();
W = W.lumpM();

% functions
kappa = @(x,y) 1*~W.isInResonator([x,y]) + 0.2*W.isInResonator([x,y]);
rho_r = @(t) 0.1/(1+0.2*cos(pi/2*t));
rho = @(x,y,t) 1*~W.isInResonator([x,y]) + rho_t(t)*W.isInResonator([x,y]);
kappa1 = 0.2;
g = @(x,y,t) sin(omega*(x - t));
u0 = @(x,y) x*0;
u1 = u0;

% matrices 
[A, M] = W.getGlobMat([1, 1/rho_r(0)], [1, kappa1]);
