%op_CSICombineCoils.m
%
% Combines the coils of MRSI data. First, coil's phases are adjusted to
% match the phase of the first coil. Then, sensetivity weights for coils are
% calculated by the formula w(i) = S(i)/sqrt(sum(S^2)) where S is intensity of the
% first point of the fids in the center of k space and i is the coil
% number. Coils are then summed by F = sum(W(i)*f(i)) where W(i) is coil's
% sensetivity weight at coil i, f(i) is the fids at coil i and F is the
% final summed signal.
%
% USAGE:
% in = op_CSICombineCoils(in)
%
% INPUT:
% in            = input MRSI object for coils to be combined
% samplePoint(optional)   = samplePoint from fid to be used for phase
%               correction. Default to first point.
% phaseMap(optional)   = 3D matrix of phases for all coils and points. Dimensions are
%                       [num coils, num x, num y]
% weightMap(optional)   = 3D matrix of weights for all coils and points.
%                         Dimensions are the same as above.
%
% OUTPUT:
% in            = MRSI object with combined coils in fids and specs.

function [MRSIStruct, phaseMap, weightMap] = op_CSICombineCoils(MRSIStruct, ...
                                        samplePoint, phaseMap, weightMap)
    arguments
        MRSIStruct (1,1) struct
        samplePoint (1,1) double = 1
        phaseMap double = []
        weightMap double = []
    end
    checkArguments(MRSIStruct);

    [MRSIStruct, pastPermute, pastSize] = reshapeDimensions(MRSIStruct, {'t', 'coils', 'y', 'x'});
    data = getData(MRSIStruct);
    if isempty(phaseMap)
        %phase map arranged by coils, y coordinate, x coordinate
        
        %calculate the angle of each each coordinate for all the coils
        phaseMap = squeeze(angle(data(samplePoint, :, :, :, :)));
    end
    %apply phase map to data
    data = applyMap(MRSIStruct, exp(-1i*phaseMap), data);

    if isempty(weightMap)
        %get weights from all positions and coils
        weights = squeeze(data(samplePoint, :, :, :, :));
        %Get the root sum squared value along coil dimension
        weightMap = normalize(weights, 1, 'norm', 2); 
    end
    %adding weights to each coil
    data = applyMap(MRSIStruct, weightMap, data);
    %add coils together
    
    dimCoils = getDimension(MRSIStruct, 'coils');
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, pastPermute, pastSize);
    data = getData(MRSIStruct);
    data = squeeze(sum(data, dimCoils));
    MRSIStruct = setData(MRSIStruct, data);

    MRSIStruct = removeDimension(MRSIStruct, 'coils');

    MRSIStruct.flags.addedrcvrs = 1;
end


%check arguments to make sure function runs as expected
function checkArguments(in)
    %some pre condition checks and setting default values
    if in.flags.addedrcvrs == 1
        error('coils already combined!')
    end
    if getFlags(in, 'spatialFT') ~= 1
        error('Please us op_CSIFourierTransform along the spatial dimension')
    end
    if getFlags(in, 'spectralFT') ~= 0
        error('Please keep the spectrum in the time domain')
    end
end



%apply map to data
function data = applyMap(MRSIStruct, map, data)
    %mutliply phases
    tSize = getSizeFromDimensions(MRSIStruct, {'t'});
    mapWithTime = repmat(map, [ones(1, ndims(map)) tSize]);

    dimensions = 1:ndims(mapWithTime);
    timeFirstOrder = circshift(dimensions, 1);
    timeFirstMap = permute(mapWithTime, timeFirstOrder);  
    data = data.*timeFirstMap;
end
