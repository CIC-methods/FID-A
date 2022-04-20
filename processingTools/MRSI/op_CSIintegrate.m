%op_CSIIntegrate.m
% Brenden Kadota, Sunnybrook 2021.
%
% USAGE:
% in=op_CSIintegrate(in);
%
% DESCRIPTION:
% Integrates the spectrum from ppmmin to ppmmax in all csi voxels.
%
% INPUTS:
% MRSIStruct        = CSI FID-A data structure
% pmmmin    = lower integration bounds
% ppmmax    = upper integration bounds
% mode      = mode (optional):
%                        -'re' (integral performed on real part (default)).
%                        -'im' (integral performed on imaginary part).
%                        -'mag' (integral performed on magnitude part).
%
% OUTPUTS:
% map       = map of integrated area
function map = op_CSIintegrate(MRSIStruct, ppmmin, ppmmax, mode)
    arguments
        MRSIStruct (1, 1) struct
        ppmmin (1, 1) double
        ppmmax (1, 1) double
        mode (1, :) char {mustBeMember(mode, {'re', 'im', 'mag'})} = 're'
    end
    % check MRSI Struct
    checkArguments(MRSIStruct);
    
    % resahpe to time, y and x dimensions
    MRSIStruct = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
    % intalize map size
    map = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));

    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'y'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'x'})
                voxel = op_CSItoMRS(MRSIStruct, x, y, "Extra", e);
                map(y, x, e) = op_integrate(voxel, ppmmin, ppmmax, mode);
            end
        end
    end
end

% argument checks. Check if spectral and spatial fourier transform have been done.
function checkArguments(in)
    if(in.flags.spectralFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spectral dimension');
    end
    if(in.flags.spatialFT == 0)
        error('FID-A Error: Input type invalid. Please fourier transform along the spacial dimension');
    end
end
