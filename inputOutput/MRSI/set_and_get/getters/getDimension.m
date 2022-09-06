function dimNumber = getDimension(MRSIStruct, dimLabel)
    arguments
        MRSIStruct (1,1) struct
        dimLabel char
    end
    switch(lower(dimLabel))
        case 'x'
            dimNumber = MRSIStruct.dims.x;
        case 'y'
            dimNumber = MRSIStruct.dims.y;
        case 'z'
            dimNumber = MRSIStruct.dims.z;
        case 'coils'
            dimNumber = MRSIStruct.dims.coils;
        case 'averages'
            dimNumber = MRSIStruct.dims.averages;
        case 'average'
            dimNumber = MRSIStruct.dims.averages;
        case 'kx'
            dimNumber = MRSIStruct.dims.kx;
        case 'ky'
            dimNumber = MRSIStruct.dims.ky;
        case 'kz'
            dimNumber = MRSIStruct.dims.kz;
        case 'subspecs'
            dimNumber = MRSIStruct.dims.subSpecs;
        case 'subspec'
            dimNumber = MRSIStruct.dims.subSpecs;
        case 'extras'
            dimNumber = MRSIStruct.dims.extras;
        case 'extra'
            dimNumber = MRSIStruct.dims.extras;
        case 't'
            dimNumber = MRSIStruct.dims.t;
        otherwise
            error('Dimension not found %s', dimLabel);
    end
end
