% author: Mauro Morini
% last modified: 06.06.24
classdef WaveGuide
    % class of rectangular 2d wave guide mesh with properties and methods to create
    % different metamaterials with resonators 
    %
    % Properties:
    % bbox: (2,2) matrix describing the rectangle domain [x1,x2;y1,y2]
    %       where x is the lower left bounding point and y the upper right
    % H: scalar number denoting the meshsize in the background (also
    %         the biggest meshsize)
    %
    % h: scalar number denoting the meshsize of the resonators
    % 
    % gd: (s, nFace) matrix each column denotes a face which add up to the
    %       total mesh, the first column always denotes the background face
    %       s is the maximal length all the face column vectors, the
    %       remaining faces with less entries are stocked up with zeros
    %       (for column description check mathworks gs decsg)
    %
    % sf: string set formula, (here just a sum of all)
    %
    % ns: (nFace,1) string vector with string denotations for Faces
    %
    % geom: geometry matrix resulting from decsg(gd,sf,ns)
    %
    % model: pde model object created by generateMesh 
    % p : (2,nP) coordinate matrix with points in columns
    % t : (4,nE) connectivity matrix with elements in rows last row is
    %       index of face to which element berlongs
    %      
    % e : (7,nE) edge matrix containing two connecting edge point indices in
    %       each column (first two rows), of the boundary edges of the domain
    %
    % Ab: (nP,nP) stiffness matrix of background 
    % Ares: (nP,nP) stiffness matrix of resonators 
    % Mb: (nP,nP) mass matrix of background 
    % Mres: (nP,nP) mass matrix of resonators 
    % MisLumped: logical, true if M is lumped

    properties  
        % geometry
        bbox 
        gd 
        sf 
        ns
        geom

        % mesh
        H = 0.5;
        h = 0.1;
        hGrad = 1.5;
        model
        p
        e
        t

        % matrices
        Ab
        Ares
        Mb
        Mres
        MisLumped
    end

    methods
        % Constructor
        function obj = WaveGuide(varargin)
            if nargin == 0
                % default constructer, no resonator
                obj.bbox = [0, 0; 20, 10];
                obj.gd = [3, 4, 0, 20, 20, 0, 0, 0, 10, 10]';
                obj.sf = "B";
                obj.ns = ["B"];
                obj = updateModel(obj);
                
            elseif nargin == 1
                % only input bounding box
                obj.bbox = varargin{1};
                obj.gd = [3, 4, obj.bbox(1,1), obj.bbox(2,1), obj.bbox(2,1), obj.bbox(1,1),...
                            obj.bbox(1,2), obj.bbox(1,2), obj.bbox(2,2), obj.bbox(2,2)]';
                obj.sf = "B";
                obj.ns = ["B"];
                obj = updateModel(obj);
            elseif nargin == 6
                % create wave guide with predefined periodic resonator subdomains
                % based on last input name
                % 
                % Inputs:
                % p0: (1,2) vector with x and y coordinate of leftmost
                %       point, startingpoint for periodicity
                % H: scalar background meshsize (max meshsize)
                % h: scalar resonator meshsize 
                % hGrad: scalar change rate between meshsizes (default 1.5)
                % [xNum,yNum]: (1,2) vector with number of subdomains in x
                %               and y direction
                % [xLen, yLen]: (1,2) vector with length of each subdomain
                %                in the respective dimension
                % r: radius of circular resonators in subdomains

                p0 = varargin{1};
                H = varargin{2}(1);
                h = varargin{2}(2);
                hGrad = varargin{2}(3);
                xNum = varargin{3}(1);
                yNum = varargin{3}(2);
                xLen = varargin{4}(1);
                yLen = varargin{4}(2);
                r = varargin{5};

                switch varargin{6}
                    case "1_in_Square"
                        obj = WaveGuide([p0; p0(1) + xNum*xLen, p0(2) + yNum*yLen]);
                        c0 = [p0(1)+xLen/2, p0(2)+yLen/2];
                        obj = setMeshsize(obj,H,h,hGrad);
                        for i = 1:xNum
                            for j = 1:yNum
                                obj = addCircRes(obj, c0+[(i-1)*xLen, (j-1)*yLen], r);
                            end
                        end
                    case  "2_in_Square"
                        obj = WaveGuide([p0; p0(1) + xNum*xLen, p0(2) + yNum*yLen]);
                        c0 = [p0(1)+xLen/2 - 1.2*r, p0(2)+yLen/2];
                        c1 = [p0(1)+xLen/2 + 1.2*r, p0(2)+yLen/2];
                        obj = setMeshsize(obj, H,h,hGrad);
                        for i = 1:xNum
                            for j = 1:yNum
                                obj = addCircRes(obj, c0+[(i-1)*xLen, (j-1)*yLen], r);
                                obj = addCircRes(obj, c1+[(i-1)*xLen, (j-1)*yLen], r);
                            end
                        end
                    case "2_in_Quadri"
                        obj = WaveGuide([p0 - [0, yNum*yLen/2]; p0 + [xNum*xLen, yNum*yLen/2]]);
                        obj = setMeshsize(obj,H,h,hGrad);
                        for i = 1:xNum
                            for j = 1:yNum
                                c = p0 + (i-1)*[xLen, yLen]/2 + (j-1)*[xLen, -yLen]/2;
                                obj = addCircRes(obj, c + [xLen/3, 0], r);
                                obj = addCircRes(obj, c + [2*xLen/3, 0], r);
                            end
                        end
                    case "honeycomb"
                        obj = WaveGuide([p0 - [0, yNum*yLen/2]; p0 + [xNum*xLen, yNum*yLen/2]]);
                        obj = setMeshsize(obj, H,h,hGrad);
                        for i = 1:xNum
                            for j = 1:yNum
                                c = p0 + (i-1)*[xLen, yLen]/2 + (j-1)*[xLen, -yLen]/2;
                                c1 = c + [xLen/3 + 3*r, 0];
                                c2 = c + [xLen/3,0] + 3*r*[cos(2*pi/3), sin(2*pi/3)];
                                c3 = c + [xLen/3,0] + 3*r*[cos(4*pi/3), sin(4*pi/3)];
                                c4 = c + 2*[xLen/3,0] + 3*r*[cos(pi/3), sin(pi/3)];
                                c5 = c + 2*[xLen/3,0] - 3*r*[xLen/3, 0];
                                c6 = c + 2*[xLen/3,0] + 3*r*[cos(5*pi/3), sin(5*pi/3)];
                                obj = addCircRes(obj, c1, r);
                                obj = addCircRes(obj, c2, r);
                                obj = addCircRes(obj, c3, r);
                                obj = addCircRes(obj, c4, r);
                                obj = addCircRes(obj, c5, r);
                                obj = addCircRes(obj, c6, r);
                            end
                        end
                    otherwise
                        error("constructed for last input: " + varargin{6} +" has not been implemented")
                end
                obj = updateModel(obj);
                
            else 
                error("Constructer for "+nargin+" inputs not implemented")
            end
        end

        % getters
        function [p,e,t] = getPet(obj)
            % returns adapted p,e,t matrices so only the actual points
            % remain (which is how they are used mostly here)
            p = obj.p';
            t = obj.t';
            t = t(:, 1:3);
            e = obj.e';
            e = e(:, 1:2);
        end

        function [A, M] = getGlobMat(obj, cA, cM)
            % returns the global stiffness and mass matrix scaled for convenience
            % Inputs:
            % cA: (2,1) matrix with scaling values for A
            % cM: (2,1) matrix with scaling values for M  

            if nargin == 0
                cA = ones(2,1);
                cM = cA;
            end
    
            A = getGlobStiffness(obj, cA);
            M = getGlobMass(obj, cM);
            end

        function A = getGlobStiffness(obj, cA)
            % returns scaled global stiffness matrix for convenience
            % Inputs:
            % cA: (2,1) matrix with scaling values for A
             if nargin == 0
                cA = ones(2,1);
             end
             A = cA(1)*obj.Ab + cA(2)*obj.Ares;
        end

        function M = getGlobMass(obj, cM)
            % returns scaled global mass matrix for convenience
            % Inputs:
            % cM: (2,1) matrix with scaling values for M
             if nargin == 0
                cM = ones(2,1);
             end
             M = cM(1)*obj.Mb + cM(2)*obj.Mres;
        end

        function nRes = getNRes(obj)
            % returns amount of resonators in mesh
            nRes = size(obj.gd, 2) - 1;
            if nRes == -1
                nRes = 0;
            end

            if nRes < 0
                error("There cannot be a negative amount of resonators, fix bug")
            end
        end
        
        function resIdx = getPointsInRes(obj)
            % index vetor with point indices of all points p which are in
            % the resonators
            resIdx = (obj.t(4,:) > 1);
            [p,e,t] = obj.getPet;
            resIdx = [t(resIdx,1);t(resIdx,2);t(resIdx,3)];
            resIdx = sort(resIdx);
            resIdx = unique(resIdx);
        end

        % setters
        function obj = setMeshsize(obj, H, h, hGrad)
            % updates meshsize parameters, depending on how many input
            % variables there are

            if nargin == 0
                % reset parameters to default
                obj.H = 0.5;
                obj.h = 0.1;
                obj.hGrad = 1.5;
            end

            if exist('H','var') 
                obj.H = H;
            end
            if exist('h','var') 
                obj.h = h;
            end
            if exist('hGrad','var') 
                obj.hGrad = hGrad;
            end
        end

        % Matrix assembly 
        function obj = assembleMatrices(obj)
            resIdx = (obj.t(4,:) > 1);
            [p,e,t] = getPet(obj);
            obj.Ab = FEM2D.stiffnessMatrix2D(p,t(~resIdx,:));
            obj.Ares = FEM2D.stiffnessMatrix2D(p,t(resIdx,:));
            obj.Mb = FEM2D.massMatrix2D(p,t(~resIdx,:));
            obj.Mres = FEM2D.massMatrix2D(p,t(resIdx,:));
            obj.MisLumped = false;
        end
        
        function obj = lumpM(obj)
            % lumps both mass matrices
            
            if obj.MisLumped
                disp("Mass matrices have already been lumped")
                return
            end

            obj.Mb = diag(sum(obj.Mb, 2));
            obj.Mres = diag(sum(obj.Mres, 2));
            obj.MisLumped = true;
        end

        % compute resonance frequency
        function eig = computeResonanceFrequencies(obj, cA, cM, nEig, sigma)
            [A,M] = getGlobMat(obj, cA, cM);
            eig = sqrt(eigs(A,M, nEig, sigma));
        end

        % adapt and update mesh
        function obj = updateModel(obj)
            % updates object and recreates mesh with current parameters
            obj.model = createpde;
            obj.geom = decsg([obj.gd],[obj.sf],[obj.ns]);
            geometryFromEdges(obj.model, obj.geom);
            generateMesh(obj.model, "GeometricOrder","linear", Hmax=obj.H,Hface={2:size(obj.gd,2), obj.h},Hgrad=obj.hGrad);
            [obj.p,obj.e,obj.t] = meshToPet(obj.model.Mesh);
        end

        function plotMesh(obj)
            % plots mesh and geometry
            figure(1)
            pdegplot(obj.geom,"FaceLabels","on")
            figure(2)
            pdeplot(obj.model)
        end
        
        % Add Resonators
        function obj = addCircRes(obj, cent, r)
        % adds a circular resonator to the grid the resonator has to be
        % fully included in the background material (with some space)
        % 
        % Inputs:
        % cent: (2,1) point vector of circle radius
        % r: scalar positive radius of circle (not bigger than bbox)
            [n, i] = size(obj.gd);
            res = [1, cent(1), cent(2), r]';
            res = [res; zeros(n-4,1)];
            obj.gd = [obj.gd, res];
            obj.ns = [obj.ns; "CR"+i];
            obj.sf = obj.sf + "+CR" + i;
            % obj = obj.updateModel();
        end

        function obj = addRectRes(obj, p)
        % adds a rectangular resonator to the grid the resonator has to be
        % fully included in the background material (with some space)
        % 
        % Inputs:
        % p: (2,2) matrix of bounding box (first row is lower left and
        %       second is upper right point
            [n, i] = size(obj.gd);
            res = [3,4,p(1,1),p(2,1),p(2,1),p(1,1),p(1,2),p(1,2),p(2,2), p(2,2)]';
            if n > 10
                res = [res; zeros(n-10,1)];
            elseif n < 10
                obj.gd = [obj.gd; zeros(10-n,i)];
            end
            obj.gd = [obj.gd, res];
            obj.ns = [obj.ns; "RR"+i];
            obj.sf = obj.sf + "+RR" + i;
            % obj = obj.updateModel();
        end

        function obj = addPolyRes(obj, p)
            % adds a polygonial resonator to the grid the resonator has to be
            % fully included in the background material (with some space)
            % 
            % Inputs:
            % p: (nP,2) matrix of nodes (first column are x and second are y
            %       values of the nodes
            [n, i] = size(obj.gd);
            m = size(p,1);
            res = [2;m;p(:,1);p(:,2)];
            m = length(res);
            if n > m
                res = [res; zeros(n-m,1)];
            elseif n < m
                obj.gd = [obj.gd; zeros(m-n,i)];
            end
            obj.gd = [obj.gd, res];
            obj.ns = [obj.ns; "PR"+i];
            obj.sf = obj.sf + "+PR" + i;
            obj = obj.updateModel();
        end

        % Resonator specific functions
        function z = isInResonator(obj, p)
            % checks if given 2d points are inside of the resonators
            % Inputs:
            % p: (nP,2) point matrix to be evaluated (not the same as the 
            %       point matrix of the mesh
            %
            % Outputs:
            % z: (nP,1) locical vector with true in all indices, for which
            %       the points in p are inside a resonator
            
            % Initializations
            nP = size(p, 1);
            nRes = getNRes(obj);
            z = false(nP,1);

            if nRes == 0
                return
            end

            res = obj.gd(:,2:end);
            for i = 1:nRes
                if res(1,i) == 1
                    r = res(4,i);
                    x = res(2,i);
                    y = res(3,i);
                    diffPRes = vecnorm(p - [x, y], 2, 2);
                    zLoc = diffPRes <= r; 
                else
                    nResLoc = res(2,i);
                    x = res(3:nResLoc+2,i);
                    y = res(nResLoc+3:(nResLoc*2+2),i);
                    zLoc = inpolygon(p(:,1), p(:,2), x, y);
                end
                z = z + zLoc;
            end

            if sum(z > 1) > 0
                warning("points are in multiple resonators at once")
            end
            z = z > 0;
        end
    end
end