%op_CSIgetLW
% Jamie Near, Sneha Senthil, Sunnybrook 2022
% USAGE:
% [FWHM]=op_CSIgetLW(in,Refppmmin,Refppmmax,zpfactor);
% 
% DESCRIPTION:
% CSI version of op_getLW
% 
% INPUTS:
% MRSIStruct         = input spectrum in structure format.
% Refppmmin  = Min of frequency range (ppm) in which to search for reference peak.
%                  (Optional.  Default = 4.4 ppm);
% Refppmmax  = Max of frequency range to (ppm) in which search for reference peak
%                  (Optional.  Default = 5/3 ppm per Tesla B0);
% zpfactor   = zero-padding factor (used for method 1.)
%                  (Optional.  Default = 8);
%
% OUTPUTS:
% FWHM       = Estimated linewidths of the MRSI struct (in Hz).

function [lw]=op_CSIgetLW(MRSIStruct,Refppmmin,Refppmmax,zpfactor);
arguments
    MRSIStruct (1, 1) struct
    Refppmmin (1, 1) double = 4.4
    Refppmmax (1, 1) double = 5.0
    zpfactor (1, 1) double =2
end

checkArguments(MRSIStruct);
[MRSIStruct, prevPermute, prevShape] = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
 lw = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));

    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                mrs = op_CSItoMRS(MRSIStruct, x, y, 'extraIndex', e); 
                LW = op_getLW(mrs, Refppmmin, Refppmmax, zpfactor,true);  %true means suppress plots.           
                lw(y, x, e) = LW;
              
                end
            end
        end
    
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);

    end

function checkArguments(MRSIStruct)
    if(getFlags(MRSIStruct, 'spatialFT') == 0)
        error('Please fourier transform along the spectral dimension')
    end
end

