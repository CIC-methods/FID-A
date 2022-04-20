classdef Voxel
    properties
        vertecies
        %Edges connects points by indexes ie (1,2) connects point 1 to
        %point 2
        edges
        %center: Center point of the voxel
        center
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
                obj.voxelFaceCenter= cat(2, coronal, sagital, axial);
                obj.center = mean(obj.voxelFaceCenter, 2);
                obj.vertecies = Voxel.getVertecies(obj.center, coronal, sagital, axial);
                obj.edges = [1,2; 2,3; 3,4; 4, 1; 1,5; 2,6; 3,7; 4,8; 5,6; 6,7; 7,8; 8,5];

                obj.index = index;
                obj.minimumCoordinates = min(obj.vertecies, [], 2);
                obj.maximumCoordinates = max(obj.vertecies, [], 2);
            end
        end

        % find_intersection
        % Finds the points that are intersecting with the plane. The points
        % are then given in sequential order to be plotted, with the first
        % point beign given a second time.
        function allIntersections = find_intersection(obj, plane, planeCoordinate)
            allIntersections = zeros(3, 4);

            pointCounter = 1;
            if(strcmpi(plane, 'coronal'))
                dimension = 2;
            elseif(strcmpi(plane, 'sagital'))
                dimension = 1;
            elseif(strcmpi(plane, 'axial'))
                dimension = 3;
            else
                error('Dimension not found')
            end

            for iEdge = 1:size(obj.edges, 1)
                edgeVerteciesIndex = obj.edges(iEdge, :);
                vertexCoordinates = obj.vertecies(:, edgeVerteciesIndex);
                firstVertex = vertexCoordinates(:, 1);
                sercondVertex = vertexCoordinates(:, 2);
                if(firstVertex(dimension) == planeCoordinate)
                    if(~ismember(firstVertex', allIntersections', 'rows'))
                        allIntersections(:, pointCounter) = firstVertex;
                        pointCounter = pointCounter + 1;
                    end
                elseif(sercondVertex(dimension) == planeCoordinate)
                    if(~ismember(firstVertex', allIntersections', 'rows'))
                        allIntersections(:, pointCounter) = sercondVertex;
                        pointCounter = pointCounter + 1;
                    end
                else
                    minVertexCoordinate = min(firstVertex(dimension), sercondVertex(dimension));
                    maxVertexCoordinate = max(firstVertex(dimension), sercondVertex(dimension));
                    if(minVertexCoordinate < planeCoordinate  && planeCoordinate < maxVertexCoordinate)
                        edgeVector = firstVertex - sercondVertex;
                        length = (planeCoordinate - sercondVertex(dimension))/ edgeVector(dimension);
                        intersectCoordinate = edgeVector*length + sercondVertex;
                        allIntersections(:, pointCounter) = intersectCoordinate;
                        pointCounter = pointCounter + 1;
                    end
                end
            end
            
            allIntersections(:, all(allIntersections == 0)) = [];

            % if no points intersect, return empty array
            if(pointCounter == 1)
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

        % getVoxelMinimumCoordinate
        % takes an object and direction and gives the smallest coordinate in the
        % x, y, or z direction of the voxel.
        function minCoord = getVoxelMinimumCoordinate(obj, direction)
            switch(direction)
                case 'sagital'
                    minCoord = obj.minimumCoordinates(1);
                case 'coronal'
                    minCoord = obj.minimumCoordinates(2);
                case 'axial'
                    minCoord = obj.minimumCoordinates(3);
            end
        end

        function minCoord = getVoxelMaximumCoordinate(obj, direction)
            switch(direction)
                case 'sagital'
                    minCoord = obj.maximumCoordinates(1);
                case 'coronal'
                    minCoord = obj.maximumCoordinates(2);
                case 'axial'
                    minCoord = obj.maximumCoordinates(3);
            end
        end
    end
    methods(Static, Access=private)
        function vertecies = getVertecies(center, coronalCenter, sagitalCenter, axialCenter)            % Vectors pointing from center to the center of each plane
            vetorToCoronal = coronalCenter - center;
            vectorToSagital = sagitalCenter - center;
            vectorToAxial = axialCenter - center;
            % get the vertecies of one side of the voxel
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
    end
end
