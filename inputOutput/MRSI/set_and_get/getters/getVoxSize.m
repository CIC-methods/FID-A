function voxSize = getVoxSize(MRSIStruct, dimension)
    switch(lower(dimension))
        case 'x'
            voxSize = MRSIStruct.voxelSize.x;
        case 'y'
            voxSize = MRSIStruct.voxelSize.y;
        case 'z'
            voxSize = MRSIStruct.voxelSize.z;
        otherwise
            error('voxSize not found for dimension %s', dimension);
    end
end
