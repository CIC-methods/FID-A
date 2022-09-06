classdef Voxel
    properties
        vertecies
        %Edges connects points by indexes ie (1, 2) connects point 1 to
        %point 2
        edges
        %center: Center point of the voxel
        voxelCenterCoordiante
        %Index: [3,1] vector mapping back to position in MRSI object's FIDS
        %or SPECS properties
        index
        %minimum values along each x,y,z coordinate.(for speedup)
        minimumCoordinates
        %maximum values along each x,y,z coordinate. (for speedup)
        maximumCoordinates
        voxelFaceCenter

    end
    methods
        function obj = Voxel(coronal, sagital, axial, index)
            if(nargin == 0)
                % used for initalization
            else
                % center coordinate of each voxels face.
                obj.voxelFaceCenter= cat(2, coronal, sagital, axial);

                % center of the voxel 
                obj.voxelCenterCoordiante = mean(obj.voxelFaceCenter, 2);

                %vertecies in an array
                obj.vertecies = Voxel.getVerteciesFromVoxelCenters(obj.voxelCenterCoordiante, coronal, sagital, axial);
                
                % connectivity array connecting vertecies together. Rows are
                % edges and columns are indecies that connect.
                obj.edges = [1,2; 2,3; 3,4; 4, 1; 1,5; 2,6; 3,7; 4,8; 5,6; 6,7; 7,8; 8,5];
                
                % voxel index in FID-A structure
                obj.index = index;

                % minimum and maximum
                obj.minimumCoordinates = min(obj.vertecies, [], 2);
                obj.maximumCoordinates = max(obj.vertecies, [], 2);
            end
        end

        % find_intersection
        % Finds the points that are intersecting with the plane. The points
        % are then given in sequential order to be plotted, with the first
        % point beign given a second time.
        function allIntersections = find_intersection(obj, dimensionLabel, planeCoordinate)
            %initalize intersection array, first dimension is coordinate
            %positiosn (x, y, z) and second is for points that intersect.
            allIntersections = zeros(3, 4);

            % counter
            interctIndex = 1;
            % get the coorisponding dimension label for dimension
            dimension = Voxel.getDimensionNumber(dimensionLabel);

            % loop through edges
            for iEdge = 1:size(obj.edges, 1)
                % get vertecies from edge
                vertexCoordinates = obj.getEdgeVertecies(iEdge);
                firstVertex = vertexCoordinates(:, 1);
                sercondVertex = vertexCoordinates(:, 2);

                % if the first vertex lies right on the plane
                if(firstVertex(dimension) == planeCoordinate)
                    % check if vertex is already added
                    if(~ismember(firstVertex', allIntersections', 'rows'))
                        % add if not in allIntersection array
                        allIntersections(:, interctIndex) = firstVertex;
                        interctIndex = interctIndex + 1;
                    end
                % if second vertex lies right on plane
                elseif(sercondVertex(dimension) == planeCoordinate)
                    % check if vertex is already added
                    if(~ismember(sercondVertex', allIntersections', 'rows'))
                        allIntersections(:, interctIndex) = sercondVertex;
                        interctIndex = interctIndex + 1;
                    end
                % check if planecoordinate lies between first and second vertex
                else
                    % get smallest and largest coordinate along plane dimension
                    minVertexCoordinate = min(firstVertex(dimension), sercondVertex(dimension));
                    maxVertexCoordinate = max(firstVertex(dimension), sercondVertex(dimension));
                    
                    % if it does intersect
                    if(minVertexCoordinate < planeCoordinate  && planeCoordinate < maxVertexCoordinate)
                        
                        % find place of intersection along edge
                        edgeVector = firstVertex - sercondVertex;
                        length = (planeCoordinate - sercondVertex(dimension))/ edgeVector(dimension);
                        intersectCoordinate = edgeVector*length + sercondVertex;
                        allIntersections(:, interctIndex) = intersectCoordinate;
                        interctIndex = interctIndex + 1;
                    end
                end
            end

            % remove empty intersection points
            allIntersections(:, all(allIntersections == 0)) = [];

            % if no points intersect, return empty array
            if(interctIndex == 1)
                allIntersections = [];
            % if points intersect order so sequentiall form a box
            else
                allDimensions = [1, 2, 3];
                plottingDimensions = setdiff(allDimensions, dimension);

                boxCoordinates = allIntersections(plottingDimensions, :);
                centerOfBox = mean(boxCoordinates, 2);
                vectorsToVertex = boxCoordinates - centerOfBox;
                theta = atan2(vectorsToVertex(2, :), vectorsToVertex(1, :)); % angle above x axis
                [~, idx] = sort(theta);
                allIntersections = allIntersections(:, idx);
                allIntersections = [allIntersections allIntersections(:,1)];
            end
        end

        % find intersection between voxel and plane. Intersection is defined as
        % plane coordinate inbetween voxels lowest and highest coordinates.
        function isIntersect = intersect(obj, planeCoordinate, plane)
            %get dimension number use for coordinates
            dimensionNumber = Voxel.getDimensionNumber(plane);

            %get the min and max coordinate of the voxel
            lowerCoordinate = obj.minimumCoordinates(dimensionNumber);
            upperCoordinate = obj.maximumCoordinates(dimensionNumber);

            % see if it intersects with the planes coordinates
            isIntersect =  upperCoordinate >= planeCoordinate && lowerCoordinate <= planeCoordinate;
        end

        % getVoxelMinimumCoordinate
        % returns the minimum coordinate based on a plane (ie. sagital)
        function minCoordinate = getVoxelMinimumCoordinate(obj, planeLabel)
            dimension = Voxel.getDimensionNumber(planeLabel);
            minCoordinate = obj.minimumCoordinates(dimension);
        end

        function maxCoordinate = getVoxelMaximumCoordinate(obj, planeLabel)
            dimension = Voxel.getDimensionNumber(planeLabel);
            maxCoordinate = obj.maximumCoordinates(dimension);
        end

        function center = getVoxelCenter(obj)
            center = obj.voxelCenterCoordiante;
        end
    end

    methods(Access=private)
        % get the vertecies from the edge index
        function vertexCoordinates = getEdgeVertecies(obj, edgeIndex)
            % getting indecies that connect
            edgeVerteciesIndex = obj.edges(edgeIndex, :);
            % gettign vertecies from those connecting indecies
            vertexCoordinates = obj.vertecies(:, edgeVerteciesIndex);
        end
    end

    methods(Static, Access=private)

        % get the vertex positions given the center coordiante of each voxel
        % face. (ie. coronal center coordiantes coorispond to the center
        function vertecies = getVerteciesFromVoxelCenters(center, coronalCenter, ...
                sagitalCenter, axialCenter)
            % get vectors pointing in each direction
            vetorToCoronal = coronalCenter - center;
            vectorToSagital = sagitalCenter - center;
            vectorToAxial = axialCenter - center;

            % get the vertecies of one side of the voxel.
            topLeftVertex = vetorToCoronal(:, 1) + vectorToSagital(:, 1) + vectorToAxial(:, 1);
            topRightVertex = vetorToCoronal(:, 1) + vectorToSagital(:, 2) + vectorToAxial(:, 1);
            bottomLeftVertex = vetorToCoronal(:, 1) + vectorToSagital(:, 1) + vectorToAxial(:, 2);
            bottomRightVertex = vetorToCoronal(:, 1) + vectorToSagital(:, 2) + vectorToAxial(:, 2);

            % create vertecies into an array
            verteciesOfCube = cat(2, topLeftVertex, topRightVertex, bottomRightVertex, bottomLeftVertex);
            % get vertecies of other side by adding vector pointint to the other
            % side
            verteciesOfSecondCube = verteciesOfCube + 2*vetorToCoronal(:, 2);
            vertecies = cat(2, verteciesOfCube, verteciesOfSecondCube);
            % re-center vertecies
            vertecies = vertecies + center;
        end

        % takes the dimension label and give the corresponding dimension number
        % used in coordiantes
        function dimension = getDimensionNumber(dimensionLabel)
            if(strcmpi(dimensionLabel, 'coronal'))
                dimension = 2;
            elseif(strcmpi(dimensionLabel, 'sagital'))
                dimension = 1;
            elseif(strcmpi(dimensionLabel, 'axial'))
                dimension = 3;
            else
                error('Dimension not found')
            end
        end
        
    end
end
