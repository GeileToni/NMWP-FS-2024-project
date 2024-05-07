% author: Mauro Morini
% last modified 09.11.23
function L = loadVector2D(p, t, f)
% calculates global load vector for linear FEM in 2D
% Inputs : 
% p : nPx2 coordinate matrix with points in rows
% t : nEx3 connectivity matrix with elements in rows
% f : function handle 
% Outputs : 
% L : global stiffness matrix nPx1

nP = size(p, 1);
nE = size(t, 1);
L = zeros(nP, 1);

% iterate over elements
for k = 1:nE
    % Element and points
    K = t(k, :);
    p0 = p(K(1), :).';
    p1 = p(K(2), :).';
    p2 = p(K(3), :).';

    % calculate J_K elementwise jacobian of the transformation
    Jk = [p1-p0, p2-p0];
    KArea = abs(det(Jk))/2;

    % approximate f by value in the center of the element K
    pCent = mean([p0, p1, p2], 2);
    fk = f(pCent(1), pCent(2));
    % fk = mean([f(p0(1), p0(2)), f(p1(1), p1(2)), f(p2(1), p2(2))], 2);

    % Assemble 
    L(t(k,:)) = L(t(k,:)) + fk*KArea/3;
end

end