%op_CSIAverage.m
%Brenden Kadota, Sunnybrook, 2021
%
% USAGE:
% [out]=op_CSIAverage(in);
%
% DESCRIPTION:
% Averages CSI FID-A structure along average dimension.
%
% INPUTS:
% in   = Input CSI FID-A structure
%
% OUTPUTS:
% out = Averaged CSI FID-A structure


function MRSIStruct = op_CSIaddNoise(MRSIStruct, sdnoise, noise)
    arguments
        MRSIStruct (1, 1) struct
        sdnoise (1, 1) double
        noise double = setDefaultNoise(MRSIStruct, sdnoise)
    end
    data = getData(MRSIStruct);
    data = data + noise;
    MRSIStruct = setData(MRSIStruct, data);
end

function noise = setDefaultNoise(MRSIStruct, sdnoise)
    dataSize = getSize(MRSIStruct);
    noiseReal = sdnoise * randn(dataSize);
    noiseImaginary = sdnoise * randn(dataSize);
    noise = noiseReal + (1i * noiseImaginary);
end