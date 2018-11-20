classdef Mesh < handle
    
    % Documentation
    % here we read a mesh written in the Tetgen or Triangle input
    % format; the user needs to provide the full path to the prefix of
    % the mesh files (filename without .ele or .node)
    
    properties %(GetAccess='private', SetAccess='private')
        dimension
        numberOfElements
        numberOfVertices
        numberOfPoints
        dimens
        N
        Areas
        ElementList
        ElementTags
        PointList
        PointMarkerList
        BoundaryFacet
        FacetList
        ElementParameters
        PointParameters
        type
    end
    


    methods
        
        
        function obj = Mesh()
            
            obj.init();
            
        end
        
        
        function init(obj)
            obj.dimension          = 0;
            obj.numberOfElements   = 0;
            obj.numberOfVertices   = 0;
            obj.numberOfPoints     = 0;
            obj.Areas              = [];
            obj.N                  = [];
            obj.dimens             = [];
            obj.ElementList        = [];
            obj.ElementTags        = [];
            obj.PointList          = [];
            obj.PointMarkerList    = [];
            obj.BoundaryFacet      = [];
            obj.FacetList          = [];
            obj.ElementParameters  = [];
            obj.PointParameters    = [];
            obj.type               = 0;
        end
        
        
        
        
        
        function readTTMesh(obj, prefix)
            
            readElements(obj, prefix)
            readNodes(obj, prefix)
            
            if (obj.dimension == 2)
                if (obj.numberOfVertices == 4)
                    obj.type = 9; % VTK_QUAD
                elseif (obj.numberOfVertices == 9)
                    obj.type = 28; % VTK_BIQUAD
                elseif (obj.numberOfVertices == 3)
                    obj.type = 5; % VTK_TRIANGLE
                end
            elseif (obj.dimension == 3)
                if (obj.numberOfVertices == 8)
                    obj.type = 12; % VTK_HEXAHEDRON
                elseif (obj.numberOfVertices == 27)
                    obj.type = 29; % VTK_TRIQUADRA_HEXAHEDRON
                elseif (obj.numberOfVertices == 10)
                    obj.type = 24; % VTK_QUADRATIC_TETRA
                elseif (obj.numberOfVertices == 4)
                    obj.type = 10; % VTK_TETRAHEDRON
                end
            end
            
            
        end
        
        
        function makePartitioningOf(obj, prefixin, prefixout, nparts)
  
        obj.readTTMesh(prefixin);
        obj.partition(nparts);
       
        obj.ElementParameters = zeros(obj.numberOfElements, 1);
        obj.writeTTMesh(prefixout);

        end
        
        
        function partition(obj, nparts)

         assembly_time_start = tic;

         elements = obj.ElementList;
        
         nv   = obj.numberOfVertices;

            n     = obj.numberOfPoints;
            nel   = obj.numberOfElements;
            
            i=1:nv;
            j=ones(1,nv);

            % ig = [1..nv, 1..nv, ... 1..nv] nv times
            ig=repmat(i,1,nv);

            % jg = [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv] each number
            % repeated nv times
            jg=repmat(1:nv,nv,1);

            %  (:) produces a column and we need a row of indices
            jg=jg(:)';


            % [1 1 ... 1, 2 2 ... 2, ... , nv nv ... nv ] with each number
            % repeated nv^2 times, these will be the indices to choose the
            % weights for each element matrix (diffusion coefficients), place
            % them in an array and multiply the arrays with the corresponding
            % matrix
            cg=repmat(1:nel, nv^2, 1);
             
            Ig    = elements(:,ig)';
            Jg    = elements(:,jg)';
            M     = ones(length(Ig(:)), 1);    
            A     = sparse(Ig(:)+1, Jg(:)+1, M, n, n);
            
            clear i; 
            clear j; 

            clear ig; 
            clear jg; 
            clear cg;

            clear C;

            assembly_time = toc(assembly_time_start);
            fprintf('ASSEMBLY LHS TIME: %10.2f\n', assembly_time);

            %TODO:
            % You have the Laplacian matrix A already assembled from the FEM
            % mesh. Call the Metis recurcive graph partitioning in order to
            % partition the FEM mesh into nparts partitions. Store the
            % partitioning vector into obj.PointParameters property. This
            % is used during creating the VTK file so that you can later
            % visualize the partitions in Visit or Paraview.
            
            partition_time_start = tic;
            %obj.PointParameters = 
            partition_time = toc(partition_time_start);
            fprintf('PARTITIONING TIME: %10.2f\n', partition_time);
            
            %TODO:
            % Compute edgecut and size of the partitions to verify they are
            % balanced. 
            %necut = 
            %part1 = 
            %part2 = 
            %partNPARTS = 
            
            %TODO:
            % Compare Metis with multilevel method from the "KaHIP_NE"
            % package. Compare time, edgecut and partition balance
        end 
        


        
        function writeTTMesh(obj, prefix)
            %writeElements(obj, prefix)
            %writeNodes(obj, prefix)

            writeMeshToVTKFile(prefix,               ... 
                               obj.ElementList,      ...
                               obj.PointList,        ...
                               obj.ElementParameters,...
                               obj.PointParameters,  ...
                               obj.type);
        end
        
        
             



        function makeParameters(obj)
            fprintf('number of points is: %d\n', obj.numberOfPoints);
            
            obj.PointParameters(1,:) = ones(1,obj.numberOfPoints);
            
            for e=1:obj.numberOfElements
                I   = obj.ElementList(e,:);
                x   = obj.PointList(I,:);
                x_c = sum(x)/8;
                x_c2= x_c(2);
                Z   = x_c(3);
                
                if (0.4+0.1*x_c2 >= Z)
                    obj.PointParameters(1,e) = 2.0;
                elseif (0.8 - 0.2*x_c2 >= Z)
                    obj.PointParameters(1,e) = 1.5;
                else
                    obj.PointParameters(1,e) = 3.0;
                end
            end
           
        end
        
   
        
        
        
        function makeRectilinearGrid(obj, N, D)
            
            obj.init();
            
            % N = [nx; ny; nz]
            % D = [dx; dy; dz]
            
            
            Nx            = N(1);
            Ny            = N(2);
            Nz            = N(3);
            m             = N(4);

            nx            = Nx + 2*m;
            ny            = Ny + 2*m;
            nz            = Nz + 2*m;
            
            dx            = D(1);
            dy            = D(2);
            dz            = D(3);
            
            obj.N         = N;
            obj.dimens    = D;
            
            x             = 0:dx:nx*dx;
            y             = 0:dy:ny*dy;
            z             = 0:dz:nz*dz;
            
            fprintf('mesh::makeRectilinearGrid constructing mesh::PointList ...\n');
            for k=1:nz+1
                pz = (k-1)*dz;
                
                for j=1:ny+1
                    py = (j-1)*dy;
                    
                    for i=1:nx+1
                        px = (i-1)*dx;
                        
                        p = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) + i;
                        obj.PointList(p,:) = [px py pz];
                        
                    end
                    
                end
                
            end
            
            Facets= [ 0 3 2 1; 4 5 6 7; 0 4 7 3; 1 2 6 5; 0 1 5 4; 2 3 7 6] + 1;
            bfacet = 0;
            fprintf('mesh::makeRectilinearGrid constructing mesh::ElementList ...\n');
            
            nel             = nx*ny*nz;
            obj.ElementTags = zeros(nel, 3);
            obj.PointParameters  = zeros(2,(nx+1)*(ny+1)*(nz+1));

            for k=1:nz
                for j=1:ny
                    for i=1:nx
                        e                = (k-1)*nx*ny + (j-1)*nx + i;
                        v1               = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) + i;
                        v2               = (k-1)*(nx+1)*(ny+1) + (j-1)*(nx+1) + i+1;
                        v3               = (k-1)*(nx+1)*(ny+1) +  j   *(nx+1) + i+1;
                        v4               = (k-1)*(nx+1)*(ny+1) +  j   *(nx+1) + i;
                        v5               =  k   *(nx+1)*(ny+1) + (j-1)*(nx+1) + i;
                        v6               =  k   *(nx+1)*(ny+1) + (j-1)*(nx+1) + i+1;
                        v7               =  k   *(nx+1)*(ny+1) +  j   *(nx+1) + i+1;
                        v8               =  k   *(nx+1)*(ny+1) +  j   *(nx+1) + i;

                        I                = [v1 v2 v3 v4 v5 v6 v7 v8];

                        obj.ElementList(e,:) = I; %[v1 v2 v3 v4 v5 v6 v7 v8];


                        tag = [ 0 0 0 ];
                        if (i <= m)
                          tag(1) = 1; 
                        elseif (i > Nx+m)
                          tag(1) = 1;
                        end

                        if (j <= m)
                          tag(2) = 1; 
                        elseif (j > Ny+m)
                          tag(2) = 1;
                        end

                        if (k > Nz+m)
                          tag(3) = 1;
                        end

                        obj.ElementTags(e,:) = tag;

                        obj.PointParameters(2,e)  = sum(tag);

                        if (k == 1)
                           bfacet = bfacet + 1;
                           obj.BoundaryFacet(bfacet, :) = [e 1];
                           %obj.BoundaryFacet(bfacet, :) = [e I(Facets(1,:))];
                        elseif (k == nz)
                           bfacet = bfacet + 1;
                           obj.BoundaryFacet(bfacet, :) = [e 2];
                           %obj.BoundaryFacet(bfacet, :) = [e I(Facets(2,:))];
                        end

                        if (i == 1)
                           bfacet = bfacet + 1;
                           obj.BoundaryFacet(bfacet, :) = [e 3];
                           %obj.BoundaryFacet(bfacet, :) = [e I(Facets(3,:))];
                        elseif (i == nx)
                           bfacet = bfacet + 1;
                           obj.BoundaryFacet(bfacet, :) = [e 4];
                           %obj.BoundaryFacet(bfacet, :) = [e I(Facets(4,:))];
                        end

                        if (j == 1)
                           bfacet = bfacet + 1;
                           obj.BoundaryFacet(bfacet, :) = [e 5];
                           %obj.BoundaryFacet(bfacet, :) = [e I(Facets(5,:))];
                        elseif (j == ny)
                           bfacet = bfacet + 1;
                           obj.BoundaryFacet(bfacet, :) = [e 6];
                           %obj.BoundaryFacet(bfacet, :) = [e I(Facets(6,:))];
                        end

                    end
                end
            end
            

            obj.Areas            = [dx*dy dx*dy dy*dz dy*dz dx*dz dx*dz];
            obj.numberOfElements = size(obj.ElementList, 1);
            obj.numberOfPoints   = size(obj.PointList, 1);
            obj.numberOfVertices = 8;
            obj.dimension        = 3;
            obj.type             = 12;
        end
        
        
        
        
        function readElements(obj, prefix)
            
            % prefix            (input): the full path to the mesh files
            %                            prefix (filename without .ele or .node)
            
            
            % 1. read the .ele file
            % %%%%%%%%%%%%%%%%%%%%%
            % this opens a file in text 't mode for read access 'r'
            filename = strcat(prefix, '.ele');
            
            fprintf('reading mesh file %s\n', filename);
            
            fid = fopen(filename, 'rt');
            
            % read the number of elements
            obj.numberOfElements = fscanf(fid, '%d', 1);
            
            % read the number of vertices per element
            obj.numberOfVertices = fscanf(fid, '%d', 1);
            

            % read whether or not we have element attributes;
            
            % ==============
            % ==> Remark <==
            % ==============
            % those we do not really need, but they must be here in case an
            % element file contains attributes, so that to make sure that the file
            % will be read correctly
            numberOfElementAttributes = fscanf(fid, '%d', 1);
            
            % now read the ElementList
            % ------------------------
            obj.ElementList = zeros(obj.numberOfElements, obj.numberOfVertices);
            for e=1:obj.numberOfElements
                
                % remove the first index which is a dummy one
                idummy = fscanf(fid, '%d', 1);
                
                for j=1:obj.numberOfVertices
                    % read the j_th vertex of the i_th element
                    vertex = fscanf(fid, '%d', 1);
                    
                    % and store it in Matlab indexing format (starting from 1)
                    obj.ElementList(e, j) = vertex;
                end
                
                % now remove the element attributes from the input file
                for j=1:numberOfElementAttributes
                    ddummy = fscanf(fid, '%d', 1);
                end
                
            end
            
            fprintf('closing mesh file %s\n', filename);
            fprintf('--------------------\n');
            
            % this closes the file
            fclose(fid);
            fprintf('numbering starts from %d\n', min(obj.ElementList(:)));
            %obj.ElementList = obj.ElementList -  min(obj.ElementList(:));
            
        end
        
        
        
        function readNodes(obj, prefix)
            
            % 2. read the .node file
            % %%%%%%%%%%%%%%%%%%%%%%
            
            % this opens a file in text 't mode for read access 'r'
            filename = strcat(prefix, '.node');
            
            fprintf('reading mesh file %s\n', filename);
            
            fid = fopen(filename, 'rt');
            % read the number of points
            obj.numberOfPoints = fscanf(fid, '%d', 1);
            
            % read the number of coordinates per point (problem obj.dimension)
            obj.dimension = fscanf(fid, '%d', 1);
            
            % read the number of point attributes
            numberOfPointAttributes = fscanf(fid, '%d', 1);
            
            % read whether or not the points have merkers
            havemarkers = fscanf(fid, '%d', 1);
            
            % ==============
            % ==> Remark <==
            % ==============
            % these we do not really need, since it they are usually zero;
            
            % now read the PointList
            % ----------------------
            obj.PointList = zeros(obj.numberOfPoints, obj.dimension);
            obj.PointMarkerList = zeros(obj.numberOfPoints, 1);
            for i=1:obj.numberOfPoints
                
                % remove the first index which is a dummy one
                idummy = fscanf(fid, '%d', 1);
                
                for j=1:obj.dimension
                    % read the j_th coordinate of the i_th point
                    coord = fscanf(fid, '%f', 1);
                    
                    % and store it in PointList
                    obj.PointList(i, j) = coord;
                end
                
                % now remove the point attributes from the input file
                for j=1:numberOfPointAttributes
                    ddummy = fscanf(fid, '%d', 1);
                end
                
                % and finally read the point markers
                if (havemarkers > 0)
                    marker = fscanf(fid, '%d', 1);
                    obj.PointMarkerList(i) = marker;
                end
                
            end
            
            fprintf('closing mesh file %s\n', filename);
            fprintf('--------------------\n');
            
            fclose(fid);
            
        end
        
        
        
        
        
        
        function writeElements(obj, prefix)
            
            % prefix            (input): the full path to the mesh files
            %                            prefix (filename without .ele or .node)
            
            
            % 1. read the .ele file
            % %%%%%%%%%%%%%%%%%%%%%
            % this opens a file in text 't mode for read access 'r'
            filename = strcat(prefix, '.ele');
            
            fprintf('writing mesh file %s\n', filename);
            
            fid = fopen(filename, 'wt');
            
            % read the number of elements
            numberOfElementAttributes=0;
            fprintf(fid, '%10d %10d %10d\n', obj.numberOfElements, obj.numberOfVertices, numberOfElementAttributes);
            
            
            % now read the ElementList
            % ------------------------
            for e=1:obj.numberOfElements
                
                % remove the first index which is a dummy one
                fprintf(fid, '%10d', e);
                
                for i=1:obj.numberOfVertices
                    fprintf(fid, '%10d', obj.ElementList(e, i));
                end
                
                fprintf(fid, '\n');
                
            end
            
            fprintf('closing mesh file %s\n', filename);
            fprintf('--------------------\n');
            
            % this closes the file
            fclose(fid);
            
        end
        
        
        
        
        
        function writeNodes(obj, prefix)
            
            % this opens a file in text 't mode for read access 'r'
            filename = strcat(prefix, '.node');
            
            fprintf('writing mesh file %s\n', filename);
            
            fid = fopen(filename, 'wt');
            % write the number of points
            fprintf(fid, '%10d', obj.numberOfPoints);
            
            % read the number of coordinates per point (problem obj.dimension)
            fprintf(fid, '%10d', obj.dimension);
            
            numberOfPointAttributes = 0;
            % read the number of point attributes
            fprintf(fid, '%10d', numberOfPointAttributes);
            
            havemarkers = 0;
            % read whether or not the points have merkers
            fprintf(fid, '%10d\n', havemarkers);
            
            for i=1:obj.numberOfPoints
                
                % remove the first index which is a dummy one
                fprintf(fid, '%10d', i);
                
                for j=1:obj.dimension
                    % read the j_th coordinate of the i_th point
                    fprintf(fid, '%23.16f', obj.PointList(i,j));
                end
                
                % and finally read the point markers
                if (havemarkers > 0)
                    fprintf(fid, '%10d', obj.PointMarkerList(i));
                end
                
                fprintf(fid, '\n');
            end
            
            fprintf('closing mesh file %s\n', filename);
            fprintf('--------------------\n');
            
            fclose(fid);
            
        end
        
        
        
        
        function neighborhood(obj)
            
            % S is the list of elements sharing each node,
            % maximum number of elements sharing a node is 8
            
            fprintf('mesh::neighborhood ...\n');
            S                    = zeros(obj.numberOfPoints, 8);
            obj.numberOfElements = size(obj.ElementList, 1);
            obj.numberOfVertices = size(obj.ElementList, 2);
            
            Index = zeros(obj.numberOfPoints, 1);
            
            for e=1:obj.numberOfElements
                
                for i=1:obj.numberOfVertices
                    V             = obj.ElementList(e, i);
                    Index(V)      = Index(V) + 1;
                    S(V,Index(V)) = e;
                end
                
            end
            
            Facets= [ 0 1 2; 0 3 1; 1 3 2; 2 3 0];
            Facets=Facets+1;
            
            
            Neighbors = zeros(obj.numberOfElements, 4); % each element may have up to 6 neighbors
            
            bindex              = 0;                        % index of the current boundary facet
            
            for e=1:obj.numberOfElements
                Vertices        = obj.ElementList(e,:)';
                for f=1:4
                    facet         = Facets(f,:);      % these are the local indices
                    
                    Facet         = Vertices(facet);  % while these are the global ones
                    
                    marker        = zeros(obj.numberOfElements, 1);
                    for v = 1:3
                        V           = Facet(v);         % for each vertex of the current facet
                        for n_i = 1:Index(V)            % loop through the surrounding elements
                            n_e       = S(V, n_i);        % and mark each element
                            marker(n_e) = marker(n_e) + 1;
                        end
                    end % v
                    
                    n_e = find(marker == 3);          % the element sharing that facet should
                    
                    %marker                           % be marked 3 times
                    if (size(n_e, 1) == 2)            % so if the number 3 has been found in
                        
                        for i=1:size(n_e,1)
                            if (n_e(i,1) ~= e)
                                n_e = n_e(i,1);
                                break;
                            end
                        end
                        
                        Neighbors(e,f) = n_e;       % the marker array the element indexed
                        
                    else                              % n_e is the f neighbor of e, otherwise
                        %fprintf('boundary facet the %d-th of element #%d\n',f,e);
                        Neighbors(e,f) = 0;             % the f-th facet is a boundary facet
                        bindex = bindex +1;
                        obj.BoundaryFacet(bindex, :) = [e f];   % the element e has its local facet f
                        obj.FacetList                = [ obj.FacetList ; Facet'];
                    end
                end % f
            end % nel
            
        end
        
        




    end % methods
    
end





