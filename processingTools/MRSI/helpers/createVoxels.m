% createVoxels.m
% creates an array of voxels based on MRSI data structure and maps back to the
% index of the data in that structure.

function vox = createVoxels(MRSIStruct)
    % get index sizes
    dimSizes = getSizeFromDimensions(MRSIStruct, {'y', 'x', 'z'});
    numberOfCoronalPoints = dimSizes(1);
    numberOfSagitalPoints = dimSizes(2);
    numberOfAxialPoints = dimSizes(3);

    %g et affine matrix
    affineMatrix = getAffineMatrix(MRSIStruct);

    % initalization
    vox(numberOfCoronalPoints, numberOfSagitalPoints, numberOfAxialPoints) = Voxel();

    for zVoxelIndex = 0:numberOfAxialPoints - 1
        for yVoxelIndex = 0:numberOfCoronalPoints - 1
            for xVoxelIndex = 0:numberOfSagitalPoints - 1
                % current voxel index
                currentIndex = [xVoxelIndex; yVoxelIndex; zVoxelIndex; 1];
                % get the six sides of the voxel
                coronal = getBoundsInWorldCoordinates('coronal', currentIndex, affineMatrix);
                sagital = getBoundsInWorldCoordinates('sagital', currentIndex, affineMatrix);
                axial = getBoundsInWorldCoordinates('axial', currentIndex, affineMatrix);
                
                % indexing by 1 in matlab. Adjust here
                currentIndex = currentIndex + 1;
                % Coronal indexing is reversed in FID-A!!! THis is to conform with
                % matlab imaging indexing for the y dimension
                currentIndex(2) = numberOfCoronalPoints - yVoxelIndex;
                vox(yVoxelIndex + 1, xVoxelIndex + 1, zVoxelIndex + 1) =...
                    Voxel(sagital, coronal, axial, currentIndex(1:3));
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
