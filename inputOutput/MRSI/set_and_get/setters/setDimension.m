function MRSIStruct = setDimension(MRSIStruct, dimLabel, value)
    switch(lower(dimLabel))
        case 'x'
            MRSIStruct.dims.x = value;
        case 'y'
            MRSIStruct.dims.y = value;
        case 'z'
            MRSIStruct.dims.z = value;
        case 'coils'
            MRSIStruct.dims.coils = value;
        case 'averages'
            MRSIStruct.dims.averages = value;
        case 'kx'
            MRSIStruct.dims.kx = value;
        case 'ky'
            MRSIStruct.dims.ky = value;
        case 'kz'
            MRSIStruct.dims.kz = value;
        case 'subspecs'
            MRSIStruct.dims.subspecs = value;
        case 'extras'
            MRSIStruct.dims.extras = value;
        case 't'
            MRSIStruct.dims.t = value;
        otherwise
            error('Dimension not found %s', dimLabel);
    end
end
