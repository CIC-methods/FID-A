% op_CSIaddPhase.m
% Brenden Kadota, Sunnybrook 2022
%
% CSI version of op_addPhase. Takes in MRSI structure, zeroth order phase,
% second order phase, and pivot to apply phase corrections. Adds phase to either
% time domain or frequency domain depending on the data dype in the data field.
%
% Input:
% MRSIStuct       = MRSI structure used in FID-A
% ph0             = First order phase in degrees (deg)
% ph1             = Second order phase in seconds (s) [default: 0]
% ppm0            = Frequency reference point (ppm) [default: 4.65]
% suppressPlot    = Plotting of frequency corrected (boolean) [default: 1]. Be
%                   careful with this option! Plots one figure for each voxel!!!
%
% Output:
% MRSIStruct      = MRSI structure used in FID-A


function MRSIStruct = op_CSIaddPhase(MRSIStruct, ph0, ph1, ppm0, suppressPlot)
arguments
    MRSIStruct (1, 1) struct
    ph0 (1, 1) double 
    ph1 (1, 1) double = 0
    ppm0 (1, 1) double = 4.65
    suppressPlot (1, 1) logical = 1
end
    checkArguments(MRSIStruct);
    [MRSIStruct, prevPermute, prevShape] = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
    data = getData(MRSIStruct);
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                mrs = op_CSItoMRS(MRSIStruct, x, y, 'extraIndex', e); 
                mrs = op_addphase(mrs, ph0, ph1, ppm0, suppressPlot);
                if(getFlags(MRSIStruct, 'spectralFT'))
                    data(:, y, x, e) = mrs.specs;
                else
                    data(:, y, x, e) = mrs.fids;
                end
            end
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);

end

function checkArguments(MRSIStruct)
    if(getFlags(MRSIStruct, 'spatialFT') == 0)
        error('Please fourier transform along the spectral dimension')
    end
end
