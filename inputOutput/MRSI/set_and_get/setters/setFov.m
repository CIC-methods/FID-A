function MRSIStruct = setFov(MRSIStruct, dimension, value)
    switch(lower(dimension))
        case('x')
            MRSIStruct.fov.x = value;
        case('y')
            MRSIStruct.fov.y = value;
        case('z')
            MRSIStruct.fov.z = value;
    end
end