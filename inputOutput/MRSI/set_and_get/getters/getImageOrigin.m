function imageOrigin = getImageOrigin(MRSIStruct, dimension)
    switch(lower(dimension))
        case 'x'
            imageOrigin = MRSIStruct.imageOrigin(1);
        case 'y'
            imageOrigin = MRSIStruct.imageOrigin(2);
        case 'z'
            imageOrigin = MRSIStruct.imageOrigin(3);
        otherwise
            error('voxSize not found for dimension %s', dimension);
    end
end