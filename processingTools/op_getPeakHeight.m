% op_getPeakHeight.m
% Jamie Near, McGill University 2017.
% 
% USAGE:
% [h]=op_getPeakHeight(in,NAAppmmin,NAAppmmax);
% 
% DESCRIPTION:
% Find the height of a peak in a spectrum.
% 
% INPUTS:
% in             = input data in matlab structure format
% NAAppmmin      = min of frequency range in which to search for peak.
%                  (Optional.  Default = 1.8 ppm (for NAA));
% NAAppmmax      = max of frequency range in which to search for peak.
%                  (Optional.  Default = 2.2 ppm (for NAA));
%
% OUTPUTS:
% h              = Peak amplitude of the desired peak.

function [h]=op_getPeakHeight(in,NAAppmmin,NAAppmmax);

%Look for NAA peak by default. 
if nargin<3
    NAAppmmax=2.2;
    if nargin<2
        NAAppmmin=1.8;
    end
end

%FIND THE NAA SIGNAL INTENSITY.  USE THE MAX PEAK HEIGHT OF THE 
%MAGNITUDE SPECTRUM INSIDE THE DESIRED SPECTRAL RANGE:
NAAwindow=in.specs(in.ppm>NAAppmmin & in.ppm<NAAppmmax);
ppmwindow=in.ppm(in.ppm>NAAppmmin & in.ppm<NAAppmmax);

maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
h=abs(NAAwindow(maxNAA_index))

