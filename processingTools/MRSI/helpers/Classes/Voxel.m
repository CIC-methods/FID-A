% Voxel.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% Voxel(position, voxelSize, rotationMatrix, fid_aIndex)
%
% DESCRIPTION:
% Class to hold voxel spatial information to plot onto mri. Could use
% refractoring to make seperate classes for edges, vertecies, and data.
%
% INPUT:    
% position          = position of voxel
% voxelSize         = size of voxel in coronal, sagital, and axial planes
% rotationMatrix    = 3x3 matrix to determine rotations
% fid_aIndex        = index to the voxel in the fid-a structure.
%
% OUTPUT:   
% Voxel             = voxel holding spatial information
classdef Voxel
    properties
        vertecies
        edges       % connectivity by using array indecies
        position
        fid_aIndex
        minimumCoordinates
        maximumCoordinates
    end
    methods
        function obj = Voxel(position, voxelSize, rotationMatrix, fidAIndex)
            if(nargin > 0)
                obj.position = position;
                obj.vertecies = Voxel.calculateVertexPositions(position, voxelSize, rotationMatrix);
                obj.edges = Voxel.getEdgesOfVoxel;
                obj.fid_aIndex = fidAIndex;
                obj.minimumCoordinates = min(obj.vertecies, [], 2);
                obj.maximumCoordinates = max(obj.vertecies, [], 2);
            end
        end

        
        function intersectionPoints = findIntersection(obj, dimensionLabel, planeCoordinate)
            %initalize intersection array
            intersectionPoints = zeros(3, 4);
            interctCounter = 1;
            % get the coorisponding dimension label for dimension
            dimension = Voxel.getDimensionNumber(dimensionLabel);

            % loop through edges
            for iEdge = 1:size(obj.edges, 1)
                % get vertecies from edge
                vertexCoordinates = obj.getVerteciesConnectedToEdge(iEdge);
                firstVertex = vertexCoordinates(:, 1);
                sercondVertex = vertexCoordinates(:, 2);

                % if the first vertex lies right on the plane
                if(firstVertex(dimension) == planeCoordinate)
                    % check if vertex is already added
                    if(~ismember(firstVertex', intersectionPoints', 'rows'))
                        % add if not in allIntersection array
                        intersectionPoints(:, interctCounter) = firstVertex;
                        interctCounter = interctCounter + 1;
                    end
                % if second vertex lies right on plane
                elseif(sercondVertex(dimension) == planeCoordinate)
                    % check if vertex is already added
                    if(~ismember(sercondVertex', intersectionPoints', 'rows'))
                        intersectionPoints(:, interctCounter) = sercondVertex;
                        interctCounter = interctCounter + 1;
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
                        intersectionPoints(:, interctCounter) = intersectCoordinate;
                        interctCounter = interctCounter + 1;
                    end
                end
            end

            % remove empty intersection points
            intersectionPoints(:, all(intersectionPoints == 0)) = [];

            % if no points intersect, return empty array
            if(interctCounter == 1)
                intersectionPoints = [];
            % if points intersect order so sequentiall form a box
            else
                allDimensions = [1, 2, 3];
                plottingDimensions = setdiff(allDimensions, dimension);

                boxCoordinates = intersectionPoints(plottingDimensions, :);
                centerOfBox = mean(boxCoordinates, 2);
                vectorsToVertex = boxCoordinates - centerOfBox;
                theta = atan2(vectorsToVertex(2, :), vectorsToVertex(1, :)); % angle above x axis
                [~, idx] = sort(theta);
                intersectionPoints = intersectionPoints(:, idx);
                intersectionPoints = [intersectionPoints intersectionPoints(:,1)];
            end
        end

        % check if there is intersection between voxel and plane.
        function isIntersect = isIntersect(obj, planeCoordinate, plane)
            %get dimension number use for coordinates
            dimensionNumber = Voxel.getDimensionNumber(plane);

            %get the min and max coordinate of the voxel
            lowerCoordinate = obj.minimumCoordinates(dimensionNumber);
            upperCoordinate = obj.maximumCoordinates(dimensionNumber);

            % see if it intersects with the planes coordinates
            isIntersect =  upperCoordinate >= planeCoordinate && lowerCoordinate <= planeCoordinate;
        end

        function minCoordinate = getVoxelMinimumCoordinate(obj, planeLabel)
            dimension = Voxel.getDimensionNumber(planeLabel);
            minCoordinate = obj.minimumCoordinates(dimension);
        end

        function maxCoordinate = getVoxelMaximumCoordinate(obj, planeLabel)
            dimension = Voxel.getDimensionNumber(planeLabel);
            maxCoordinate = obj.maximumCoordinates(dimension);
        end

        function center = get.position(obj)
            center = obj.position;
        end
    end

    % static methods
    methods(Access=private)
        % get the vertecies from the edge index
        function vertexCoordinates = getVerteciesConnectedToEdge(obj, edgeIndex)
            % getting indecies that connect
            edgeVerteciesIndex = obj.edges(edgeIndex, :);
            % gettign vertecies from those connecting indecies
            vertexCoordinates = obj.vertecies(:, edgeVerteciesIndex);
        end
    end


    methods(Static)
        % get the vertex positions given the center coordiante of each voxel
        % face. (ie. coronal center coordiantes coorispond to the center
        function vertecies = calculateVertexPositions(position, voxelSize, rotationMatrix)

            unitVoxelCoordinates = [0    0    0;
                                    1    0    0;
                                    1    1    0;
                                    0    1    0;
                                    0    0    1;
                                    1    0    1;
                                    1    1    1;
                                    0    1    1;];
            % center the unit voxel around zero
            unitVoxelCoordinates = unitVoxelCoordinates - 0.5;
            
            % multiply each by the size
            voxelCoordinatesUnrotated = voxelSize' .* unitVoxelCoordinates';
            voxelCoordinatesRotated = rotationMatrix * voxelCoordinatesUnrotated;
            vertecies = voxelCoordinatesRotated + position;
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

        function edges = getEdgesOfVoxel
            edges = [1 2
                     2 3
                     3 4
                     1 4
                     1 5
                     2 6
                     3 7
                     4 8
                     5 6
                     6 7
                     7 8
                     5 8];
        end
    end
end