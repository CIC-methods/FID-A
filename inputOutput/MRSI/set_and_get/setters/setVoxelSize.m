function MRSIStruct = setVoxelSize(MRSIStruct, dimension, value)
    switch(lower(dimension))
        case('x')
            MRSIStruct.voxelSize.x = value;
        case('y')
            MRSIStruct.voxelSize.y = value;
        case('z')
            MRSIStruct.voxelSize.z = value;
    end
end