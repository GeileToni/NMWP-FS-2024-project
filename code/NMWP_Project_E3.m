% author: Mauro Morini  
% last modified: 24.04.24
clc;clear;close all;

% Initializations
H = 2.^(-(2:6))/2;
r = 0.5;
T = [0.25, 0.5, 0.75, 1];
L2err = zeros(1, length(H));

% functions
u0 = @(x,y) 1/pi*cos(pi*x);
u1 = @(x,y) sin(pi*x);
g = @(x,y,t) sin(pi*(x - t));
uexact = @(x,y,t) 1/pi*cos(pi*(x-t));

for l = 1:length(H)
    h = H(l);
    dt = r*h;
    % mesh
    [p, t, e] = generateMesh2dUnitSquare(h);
    
    % calculate numerical solution
    [uh_T2, tvec] = waveEqLF2D(u0, u1, g, T, dt, p, e, t);

    % exact solution at tvec(end) ~= T(end)
    Uexact = uexact(p(:,1),p(:,2),tvec(end));

    % error
    L2err(l) = error2d(p,t,@(x,y) uexact(x,y,tvec(end)),uh_T2(:,end));
end

%% plots

% plot error
figure(1)
plog = loglog(H, L2err, H, H.^2, '--');
legend("L2-error", "h^2")
plog(1).LineWidth = 1;
plog(2).LineWidth = 1;
xlabel("h")
ylabel("error")


% plot exact and numerical solution in finest mesh
figure(2);
tld = tiledlayout('flow');
title(tld, "numerical solution of " + ...
    "$\frac{1}{\pi} cos(\pi(x - t))$ with FEM", 'Interpreter', 'latex')
for i = 1:length(T)
    nexttile
    trisurf(t,p(:,1),p(:,2),uh_T2(:,i), 'EdgeColor', 'none')
    %view(2)
    title("t = " + T(i))
    xlabel("x")
    ylabel("y")
    zlabel("uh")
end

figure(3);
tld = tiledlayout('flow');
title(tld, "exact solution of " + ...
    "$\frac{1}{\pi} cos(\pi(x - t))$ with FEM", 'Interpreter', 'latex')
for i = 1:length(T)
    nexttile
    trisurf(t,p(:,1),p(:,2),uexact(p(:,1),p(:,2),T(i)), 'EdgeColor', 'none')
    %view(2)
    title("t = " + T(i))
    xlabel("x")
    ylabel("y")
    zlabel("uh")
end

