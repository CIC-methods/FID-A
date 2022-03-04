function coordinates = getCoordinates(MRSIStruct, dimension)
    switch(lower(dimension))
        case 'x'
            coordinates = MRSIStruct.coordinates.x;
        case 'y'
            coordinates = MRSIStruct.coordinates.y;
        case 'z'
            coordinates = MRSIStruct.coordinates.z;
        case 'kx'
            [fovKx, deltaKx] = getKSpaceValues(MRSIStruct, 'x');
            coordinates = createCoordinates(fovKx/2, deltaKx);
        case 'ky'
            [fovKy, deltaKy] = getKSpaceValues(MRSIStruct, 'y');
            coordinates = createCoordinates(fovKy/2, deltaKy);
        case 'kz'
            [fovKz, deltaKz] = getKSpaceValues(MRSIStruct, 'z');
            coordinates = createCoordinates(fovKz/2, deltaKz);
        otherwise
            error('fov not found for dimension %s', dimension);
    end
end

function [kFov, kDelta] = getKSpaceValues(MRSIStruct, dimension)
    imageFov = getFov(MRSIStruct, dimension);
    kDelta = 1/imageFov;
    imageCoordinates = getCoordinates(MRSIStruct, dimension);
    if(length(imageCoordinates) == 1)
        kFov = 1/imageFov;
    else
        kFov = 1/(imageCoordinates(2) - imageCoordinates(1));
    end
end