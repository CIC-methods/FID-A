%op_CSIIntegrate.m
%Brenden Kadota, Sunnybrook 2021.
%
% USAGE:
% in=op_CSIintegrate(in);
%
% DESCRIPTION:
% Integrates the spectrum from ppmmin to ppmmax in all csi voxels.
%
% INPUTS:
% in        = CSI FID-A data structure
% pmmmin    = lower integration bounds
% ppmmax    = upper integration bounds
% mode      = mode (optional):
%                        -'re' (integral performed on real part (default)).
%                        -'im' (integral performed on imaginary part).
%                        -'mag' (integral performed on magnitude part).
%
% OUTPUTS:
% out       = map of integrated area
function [map, signal, noised] = op_CSIgetSNR(MRSIStruct, NAAppmmin, NAAppmmax, noiseppmmin, noiseppmmax)
    arguments
        MRSIStruct (1, 1) struct
        NAAppmmin (1, 1) double = nan
        NAAppmmax (1, 1) double = nan
        noiseppmmin (1, 1) double = nan
        noiseppmmax (1, 1) double = nan
    end
    checkArguments(MRSIStruct);
    arguments = removeNanArguments(NAAppmmin, NAAppmmax, noiseppmmin, noiseppmmax);


    MRSIStruct = reshapeDimensions(MRSIStruct, {'y', 'x', 't'});
    map = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));
    signal = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));
    noised = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                voxel = op_CSItoMRS(MRSIStruct, x, y, "Extra", e);
                [voxSNR, voxSignal, voxNoiseD] = op_getSNR(voxel, arguments{:});
                signal(y, x, e) = voxSignal;
                noised(y, x, e) = voxNoiseD;
                map(y, x, e) = voxSNR;
            end
        end
    end
end

function checkArguments(in)
    if(in.flags.spectralFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spectral dimension');
    end
    if(in.flags.spatialFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spacial dimension');
    end
end

function args = removeNanArguments(NAAppmmin, NAAppmmax, noiseppmmin, noiseppmmax)
    args = {NAAppmmin, NAAppmmax, noiseppmmin, noiseppmmax};
    emptyIndecies = zeros(1, size(args, 2), 'logical');
    for iArgument = 1:length(args)
        if(isnan(args{iArgument}))
            emptyIndecies(iArgument) = 1;
        end
    end
    args(emptyIndecies) = [];
end