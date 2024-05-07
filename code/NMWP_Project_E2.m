% author: Mauro Morini
% last modified 14.04.24
clc;clear;close all;

% Initializations
H = 2.^(-(3:6))/2;
numIter = length(H);
L2err = zeros(1, numIter);

% functions
syms x y
u = sin(pi*x)*sin(pi*y);
f = -laplacian(u) + u;
u = matlabFunction(u);
f = matlabFunction(f);
% g = @(x,y) (x == 0 || y == 0).*(-pi).*sin(x*pi).*sin(y*pi) + (x == 1 || y == 1).*pi.*sin(x*pi).*sin(y*pi);
g = @(x,y)  (x == 1).*(pi).*cos(x*pi).*sin(y*pi) + ...
            (y == 1).*(pi).*sin(x*pi).*cos(y*pi) + ...
            (x == 0).*(-pi).*cos(x*pi).*sin(y*pi) + ...
            (y == 0).*(-pi).*sin(x*pi).*cos(y*pi);

for i = 1:numIter
    [p, t, e] = generateMesh2dUnitSquare(H(i));
    uh = solve_elliptic_BVP_2d_FEM_Neumann(p,t,e,f,g);
    u_exact = u(p(:,1), p(:, 2));
    L2err(i) = error2d(p, t, u, uh);
end

%% plot
figure()
plog = loglog(H,L2err,H,H.^2,'--');
legend("L2-error", "h^2")
title("error of linear FEM")
xlabel("h")
plog(1).LineWidth = 1;
plog(2).LineWidth = 1;

% figure()
% trisurf(t,p(:,1),p(:,2),u_exact, 'EdgeColor', 'none')
figure()
trisurf(t,p(:,1),p(:,2),uh, 'EdgeColor', 'none')
xlabel("x")
ylabel("y")
zlabel("uh")
title("numerical solution")

% pdeplot(p', t')