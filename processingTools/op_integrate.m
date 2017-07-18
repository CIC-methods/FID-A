% op_integrate.m
% Jamie Near, McGill University 2016.
% 
% USAGE:
% [int]=op_integrate(in,ppmmin,ppmmax,mode);
% 
% DESCRIPTION:
% Basic peak integration over a specified frequency range.  By default,
% this function integrates under the real part of the curve, but it can
% also be made to integrate the imaginary part or the magnitude part by
% changing the "mode" parameter.
% 
% INPUTS:
% in         = input data in matlab structure format
% ppmmin     = min of frequncy range (in ppm) in which to calculate integral.
% ppmmax     = max of frequncy range (in ppm) in which to calculate integral.
% mode       = mode (optional): 
%                        -'re' (integral performed on real part (default)).
%                        -'im' (integral performed on imaginary part).
%                        -'mag' (integral performed on magnitude part).
%
% OUTPUTS:
% int        = Estimated area under the curve for the desired frequency range.


function [int]=op_integrate(in,ppmmin,ppmmax,mode);

if nargin<4
    mode='re';
    if nargin<3
        error('ERROR:  Not enough input arguments!  ABORTING!');
    end
end

%Now calculate the integral:
switch mode
    case 're'
        int=sum(real(in.specs(in.ppm>ppmmin & in.ppm<ppmmax)),1);
    case 'im'
        int=sum(imag(in.specs(in.ppm>ppmmin & in.ppm<ppmmax)),1);
    case 'mag'
        int=sum(abs(in.specs(in.ppm>ppmmin & in.ppm<ppmmax)),1);
    otherwise
        error('ERROR: Mode not recognized.  Aborting!!');
end

        
