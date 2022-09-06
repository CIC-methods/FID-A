function MRSIStruct = setCoordinates(MRSIStruct, dimension, value)
    switch(lower(dimension))
        case 'x'
            MRSIStruct.coordinates.x = value;
        case 'y'
            MRSIStruct.coordinates.y = value;
        case 'z'
            MRSIStruct.coordinates.z = value;
        otherwise
            error('fov not found for dimension %s', dimension);
    end
end
