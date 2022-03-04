function fov = getFov(MRSIStruct, dimension)
    switch(lower(dimension))
        case 'x'
            fov = MRSIStruct.fov.x;
        case 'y'
            fov = MRSIStruct.fov.y;
        case 'z'
            fov = MRSIStruct.fov.z;
        case 'kx'
            xCoordinates = getCoordinates(MRSIStruct, 'x');
            deltaX = xCoordinates(2) - xCoordinates(1);
            fov = 1/deltaX;
        case 'ky'
            yCoordinates = getCoordinates(MRSIStruct, 'y');
            deltaY = yCoordinates(2) - yCoordinates(1);
            fov = 1/deltaY;
        case 'kz'
            zCoordinates = getCoordinates(MRSIStruct, 'z');
            deltaz = zCoordinates(2) - zCoordinates(1);
            fov = 1/deltaz;
        otherwise
            error('fov not found for dimension %s', dimension);
    end
end
