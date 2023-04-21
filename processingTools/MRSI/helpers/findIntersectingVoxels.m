function [sagital_voxels, coronal_voxels, transverse_voxels] = findIntersectingVoxels(voxels, center)
    numSagitalVoxels = 1;
    numCoronalVoxels = 1;
    numTransverseVoxels = 1;
    sagital_voxels(numel(voxels)) = Voxel();
    coronal_voxels(numel(voxels)) = Voxel();
    transverse_voxels(numel(voxels)) = Voxel();
    for i = 1:numel(voxels)
        if(voxels(i).isIntersect(center(1), 'sagital'))
            sagital_voxels(numSagitalVoxels) = voxels(i);
            numSagitalVoxels = numSagitalVoxels + 1;
        end
        if(voxels(i).isIntersect(center(2), 'coronal'))
            coronal_voxels(numCoronalVoxels) = voxels(i);
            numCoronalVoxels = numCoronalVoxels + 1;
        end
        if(voxels(i).isIntersect(center(3), 'axial'))
            transverse_voxels(numTransverseVoxels) = voxels(i);
            numTransverseVoxels = numTransverseVoxels + 1;
        end
    end

    % Remove empty cells in the array
    if(numSagitalVoxels < numel(voxels))
        sagital_voxels(numSagitalVoxels:end) = [];
    end
    if(numCoronalVoxels < numel(voxels))
        coronal_voxels(numCoronalVoxels:end) = [];
    end
    if(numTransverseVoxels < numel(voxels))
        transverse_voxels(numTransverseVoxels:end) = [];
    end
end