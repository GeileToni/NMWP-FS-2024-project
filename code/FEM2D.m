% author: Mauro Morini  
% last modified: 08.05.24
classdef FEM2D
    % combination of codes for FE in 2D

    properties
        Property1
    end

    methods (Static)
        % author: Mauro Morini
        % last modified 24.04.24
        function A = stiffnessMatrix2D(p, t, c)
        % calculates global stiffness matrix for linear FEM in 2D int_omega
        % a(x)*phi_i*phi_j dx
        % Inputs : 
        % p : nPx2 coordinate matrix with points in rows
        % t : nEx3 connectivity matrix with elements in rows
        % c : if exists function handle of bilinear form
        % Outputs : 
        % A : global stiffness matrix nPxnP
        
        if ~exist('c','var') 
            c = @(x,y) 1; 
        end
        
        nP = size(p, 1);
        nE = size(t, 1);
        A = sparse(nP, nP);
        
        % iterate over elements
        for k = 1:nE
            % Element and points
            K = t(k, :);
            p0 = p(K(1), :).';
            p1 = p(K(2), :).';
            p2 = p(K(3), :).';
        
            cent = mean([p0, p1, p2], 2);       % center point
        
            % calculate J_K elementwise jacobian of the transformation
            Jk = [p1-p0, p2-p0];
        
            K_area = abs(det(Jk))/2;
            DN = [-1, -1; 1, 0; 0, 1];
        
            % elementwise local stiffness matrix
            A_loc = c(cent(1),cent(2))*K_area*DN*inv(Jk'*Jk)*DN';
        
            % Assembling
            A(t(k, :), t(k, :)) = A(t(k, :), t(k, :)) + A_loc;
        end
        end
        
        % author: Mauro Morini (adapted from "The Finite Element Theory" by
        % Larson/Bengzon chapter 4.6)
        % last modified 23.04.24
        function R = neumannMassMatrix2D(p, e, c)
        % assembles mass matrix for neumann b.c. on the edges int_gamma
        % c*phi_i*phi_j ds
        %
        % Inputs:
        % p : (nP,2) point matrix with x values in the first and y values in the
        %       second column
        % e : (nE,2) edge matrix containing two connecting edge point indices in
        %       each row, of the boundary edges of the domain
        % c : function handle from weak formulation default is 1
        %
        % Outputs:
        % R :  (nP,nP) Neumann mass matrix with R_ij = int_{e_ij} c phi_i phi_j ds
        
        if ~exist('c','var') 
            c = @(x,y) 1; 
        end
        
        % Initializations
        nP = size(p, 1);
        nE = size(e, 1);
        R = sparse(nP,nP);
        
        for i = 1:nE
            E = e(i,:);
            x1 = p(E(1),:);                               % edge point 1
            x2 = p(E(2),:);                               % edge point 2
            e_h = norm(x1 - x2, 2);                       % edge length
            cent = mean([x1', x2'], 2);
        
            % local element contributions, same as for mass matrix but here only
            % 2x2 because the integral is over e
            Rloc = e_h/6*[2, 1; 1 2]*c(cent(1), cent(2));          
            R(E, E) = R(E, E) + Rloc;
        end
        end
        
        % author: Mauro Morini (adapted from "The Finite Element Theory" by
        % Larson/Bengzon chapter 4.6)
        % last modified 16.04.24
        function L = neumannLoadVector2D(p, e, f)
        % assembles load vector for neumann b.c. : du/dn = f on boundary Gamma given by
        % the edge matrix e for 2d linear FE. In the weak formulation, the above
        % b.c. has the following impact on the LHS: int_Gamma v*f ds for a
        % testfunction v. 
        %
        % Inputs:
        % p : (nP,2) point matrix with x values in the first and y values in the
        %       second column
        % e : (nE,2) edge matrix containing two connecting edge point indices in
        %       each row, of the boundary edges of the domain
        % f : function handle of neumann b.c.
        %
        % Outputs:
        % L : (nP,1) neumann load vector containing the pointwise impact of the
        %       b.c. in each entry (zero for internal points)
        
        % Initializations
        nP = size(p, 1);
        nE = size(e, 1);
        L = zeros(nP,1);
        
        for i = 1:nE
            E = e(i,:);
            x1 = p(E(1),:);                               % edge point 1
            x2 = p(E(2),:);                               % edge point 2
            e_h = norm(x1 - x2, 2);                       % edge length
            cent = mean([x1', x2'], 2);
        
            % local element contributions
            Lloc = e_h/2*f(cent(1),cent(2))*[1;1];          
            L(E) = L(E) + Lloc;
        end
        end
    
        % author: Mauro Morini  
        % last modified: 24.04.24
        function M = massMatrix2D(p, t, c)
        % calculates global mass matrix for linear 
        % FE in 2D M(i,j) = int_Omega c(x)*phi_i*phi_j
        % Inputs : 
        % p : nPx2 coordinate matrix with points in rows
        % t : nEx3 connectivity matrix with elements in rows
        % c : if exists function handle of bilinear form
        %
        % Outputs : 
        % M : global stiffness matrix nPxnP
        
        if ~exist('c','var') 
            c = @(x,y) 1; 
        end
        
        nP = size(p, 1);
        nE = size(t, 1);
        M = sparse(nP, nP);
        
        % iterate over elements
        for k = 1:nE
            % Element and points
            K = t(k, :);
            p0 = p(K(1), :).';
            p1 = p(K(2), :).';
            p2 = p(K(3), :).';
        
            cent = mean([p0, p1, p2], 2);       % center point
        
            % calculate J_K elementwise jacobian of the transformation
            Jk = [p1-p0, p2-p0];
            K_area = abs(det(Jk))/2;
        
            % elementwise local stiffness matrix
            Mloc = [2, 1, 1; 1, 2, 1; 1, 1, 2];
            Mloc = c(cent(1),cent(2))*K_area*(1/12)*Mloc;
        
            % Assembling
            M(t(k, :), t(k, :)) = M(t(k, :), t(k, :)) + Mloc;
        end
        end
    
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
    
        % author: Mauro Morini  
        % last modified: 18.11.23
        function L2err = error2d(p, t, u, uh)
        % calculates the L^2 error of a linear FE soltion in 2D
        %
        % Inputs : 
        % p : nPx2 coordinate matrix with points in rows
        % t : nEx3 connectivity matrix with elements in rows
        % u : function handle @(x,y) of exact solution u
        % uh : nPx1 column vector of FE solution
        %
        % Output : 
        % L2err : scalar value of L^2 error (||u-uh||_{L^2(omega)})
        
        % Initalizations
        [nE, r] = size(t);
        nP = size(p, 1);
        L2err = 0;
        KhatArea = 1/2;
        
        % weights of QF
        w = (1/3)*[1, 1, 1];
        
        % nodes of QF 2xr
        xi = [0, 0; 1, 0; 0, 1]';
        
        % shape functions
        N = {@(xi) 1 - [1,1]*xi, @(xi) [1,0]*xi, @(xi) [0,1]*xi};
        
        % iterate over elements
        for i = 1:nE
            
            % Element and points
            K = t(i, :);
            p0 = p(K(1), :).';
            p1 = p(K(2), :).';
            p2 = p(K(3), :).';
        
            % calculate J_K elementwise jacobian of the transformation
            Jk = [p1-p0, p2-p0];
        
            % element map
            Fk = @(xi) p0 + Jk*xi;
        
            KArea = abs(det(Jk))/2;
        
            % calculate Quadrature
            Qval = zeros(1, r);
            for j = 1:r
                uhVal = 0;
                for q = 1:r
                    uhVal = uhVal + uh(t(i, q))*N{q}(xi(:, j));
                end
                x = Fk(xi(:, j));
                Qval(j) = abs(u(x(1), x(2)) - uhVal)^2;
            end
            L2err = L2err + (KhatArea*KArea)*Qval*w';
        end
        L2err = sqrt(L2err);
        end
    end
end