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
function map = op_CSIintegrate(MRSIStruct, ppmmin, ppmmax, mode)
    checkArguments(MRSIStruct);
    
    MRSIStruct = reshapeDimensions(MRSIStruct, {'y', 'x'});
    map = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));
    for e = getSizeFromDimensions(MRSIStruct, 'extras')
        for x = getSizeFromDimensions(MRSIStruct, {'y'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'x'})
                voxel = op_CSItoMRS(MRSIStruct, x, y, "Extra", e);
                if(exist("mode", 'var'))
                    map(y, x, e) = op_integrate(voxel, ppmmin, ppmmax, mode);
                else
                    map(y, x, e) = op_integrate(voxel, ppmmin, ppmmax);
                end

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
