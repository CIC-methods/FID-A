function MRSIStruct = op_CSIZeroFillTime(MRSIStruct, numPointsToZeroFill)
    arguments
        MRSIStruct (1, 1) struct
        numPointsToZeroFill (1, 1) double {mustBePositive}
    end
    checkArguments(MRSIStruct)
    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t'});
    data = getData(MRSIStruct);
    extraDimSize = getSizeFromDimensions(MRSIStruct, {'extras'});
    zeroFill = zeros(numPointsToZeroFill, extraDimSize);
    data = cat(1, data, zeroFill);

    for iDim = 1:size(prevPermute, 2)
        dimLabel = prevPermute{2, iDim};
        if strcmp(dimLabel, 't')
            prevSize(iDim) = prevSize(iDim) + numPointsToZeroFill;
        end
    end
    t = getSpectralTime(MRSIStruct);
    dwellTime = 1/getSpectralWidth(MRSIStruct);
    zeroFillTime = 0:dwellTime:dwellTime * (numPointsToZeroFill - 1);
    newTime = cat(2, t, zeroFillTime);
    MRSIStruct = setSpectralTime(MRSIStruct, newTime);
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);

end

function checkArguments(MRSIStruct)
    if(getFlags(MRSIStruct, 'spectralFT') == true)
        error('Please have the convert the ppm axis to time')
    end
end