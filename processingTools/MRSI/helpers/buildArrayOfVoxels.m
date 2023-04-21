% createVoxels.m
% creates an array of voxels based on MRSI data structure and maps back to the
% index of the data in that structure.

function vox = buildArrayOfVoxels(MRSIStruct)
    % get index sizes
    spatialSizes = getSizeFromDimensions(MRSIStruct, {'y', 'x', 'z'});
    numberOfCoronalPoints = spatialSizes(1);
    numberOfSagitalPoints = spatialSizes(2);
    numberOfAxialPoints = spatialSizes(3);

    affineMatrix = getAffineMatrix(MRSIStruct);
    scaleAndRotationMatrix = affineMatrix(1:3, 1:3);

    coronalVoxelSize = getVoxSize(MRSIStruct, 'y');
    sagitalVoxelSize = getVoxSize(MRSIStruct, 'x');
    axialVoxelSize = getVoxSize(MRSIStruct, 'z');
    voxelSize = [coronalVoxelSize, sagitalVoxelSize, axialVoxelSize];

    rotationMatrix = scaleAndRotationMatrix * inv(eye(3) .* voxelSize);

    % initalize 3d array holding voxels
    vox(numberOfCoronalPoints, numberOfSagitalPoints, numberOfAxialPoints) = Voxel();
    
    for zVoxelIndex = 0:numberOfAxialPoints - 1
        for yVoxelIndex = 0:numberOfCoronalPoints - 1
            for xVoxelIndex = 0:numberOfSagitalPoints - 1
                % current voxel index
                currentIndex = [xVoxelIndex; yVoxelIndex; zVoxelIndex];
                voxelPosition = affineMatrix * [currentIndex; 1];
                voxelPosition = voxelPosition(1:3);
                % get the six sides of the voxel
                
                % indexing by 1 in matlab. Adjust here
                currentIndex = currentIndex + 1;
                vox(yVoxelIndex + 1, xVoxelIndex + 1, zVoxelIndex + 1) =...
                    Voxel(voxelPosition, voxelSize, rotationMatrix, currentIndex);
            end
        end
    end
end

% get world coordinate bounds of voxel indexes
function worldCoordinates = getBoundsInWorldCoordinates(direction, ...
        currentVoxelIndex, affineMatrix)
    switch(direction)
        case('sagital')
            directionIndex = 1;
        case('coronal')
            directionIndex = 2;
        case('axial')
            directionIndex = 3;
    end
    lowerBoundIndex = currentVoxelIndex;
    upperBoundIndex = currentVoxelIndex;
    lowerBoundIndex(directionIndex) = lowerBoundIndex(directionIndex) - 0.5;
    upperBoundIndex(directionIndex) = upperBoundIndex(directionIndex) + 0.5;

    worldCoordinates = zeros(4, 2);
    worldCoordinates(:, 1) = affineMatrix * lowerBoundIndex;
    worldCoordinates(:, 2) = affineMatrix * upperBoundIndex;
    % remove fourth coordinate. Not used after affine transformation.
    worldCoordinates(4, :) = [];
end
