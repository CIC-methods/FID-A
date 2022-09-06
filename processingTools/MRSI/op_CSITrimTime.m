% op_CSIZeroFill.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% [out]=io_loadspec_rda(rda_filename);
% 
% DESCRIPTION:
% Removes the first the first time point in all fids until startIndex is reached
% 
% INPUTS:
% MRSIStruct   = FID-A CSI structure
% startIndex    = Number of points to pad MRSIStruct the x direction
%
% OUTPUTS:
% out        = Output CSI structure with zero padding

function MRSIStruct = op_CSITrimTime(MRSIStruct, startIndex)
    
    [MRSIStruct, prev_permute, prev_shape] = reshapeDimensions(MRSIStruct, {'t'});
    timeDimension = getDimension(MRSIStruct, 't');
    prev_shape(timeDimension) = prev_shape(timeDimension) - startIndex + 1; 
    data = getData(MRSIStruct);
    data = data(startIndex:end,:);
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prev_permute, prev_shape);

    newTimeLength = getSizeFromDimensions(MRSIStruct, {'t'});
    time = getAdcTime(MRSIStruct);
    newTime = time(1:newTimeLength);
    MRSIStruct = setAdcTime(MRSIStruct, newTime);
    if(getFlags(MRSIStruct, 'isCartesian') == 1)
        MRSIStruct = setSpectralTime(MRSIStruct, newTime);
    end
    if(getFlags(MRSIStruct, 'spectralft') == true)
        spectralWidth = getSpectralWidth(MRSIStruct);
        step = spectralWidth/newTimeLength;
        lb = -spectralWidth/2 + step/2;
        ub = spectralwidth/2 - step/2;

        frequency=lb:step:ub;

        gamma = 42.577478518e6;
        %calculating the ppm
        ppm = (-frequency/(MRSIStruct.Bo*gamma))*10e6;
        ppm = ppm + 4.65;

        MRSIStruct = setPPM(MRSIStruct, flip(ppm));
    end
end
