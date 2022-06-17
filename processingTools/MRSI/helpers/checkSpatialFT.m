function checkSpatialFT(MRSIStruct) 
    if(~getFlags(MRSIStruct, 'spatialFT'))
        error('Please fourier transform along the spatial dimension')
    end
end
