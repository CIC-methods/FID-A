% op_filter.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [MRSIStruct,lorenzianFilter]=op_filter(MRSIStruct,lineBroadening);
% 
% DESCRIPTION:
% Perform lMRSIStructe broadening by multiplying the time domain signal by an
% exponential decay function.  
% 
% INPUTS:
% MRSIStruct     = input data in matlab structure format.
% lineBroadening     = lMRSIStructe broadening factor in Hz.
%
% OUTPUTS:
% MRSIStruct    = Output followMRSIStructg alignment of averages.  
% lorenzianFilter    = Exponential (time domaMRSIStruct) filter envelope that was applied.

function [MRSIStruct, lorenzianFilter] = op_CSIFilterTime(MRSIStruct, lineBroadening)
    if getFlags(MRSIStruct, 'filtered')
        cont = input('WARNING: Line Broadening has already been performed!  Continue anyway?  (y or n)');
        if lower(cont) ~= 'y'
            error('STOPPING');
        end
    end

    %get t2 value from linebroading value
    t2 = 1 / (pi * lineBroadening);

    %Create an exponential decay (lorentzian filter):
    time = getTime(MRSIStruct);
    lorenzianFilter = exp(-time / t2);

    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t'});
    extraSize = getSizeFromDimensions(MRSIStruct, {'extras'});
    lorenzianFilterDataSize = repmat(lorenzianFilter, [1, extraSize]);

    %get data
    data = getData(MRSIStruct);
    %Now multiply the data by the filter array.
    data = data .* lorenzianFilterDataSize;

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeDimensions(MRSIStruct, prevPermute, prevSize);

    MRSIStruct = setFlags(MRSIStruct, 'writtentostruct', true);
    MRSIStruct = setFlags(MRSIStruct, 'filtered', true);
end
