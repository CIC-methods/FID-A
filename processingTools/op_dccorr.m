% op_dccorr.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,dcOffset]=op_dccorr(in,mode,var1);
% 
% DESCRIPTION:
% Do a DC Correction on the data.  This method is a frequency domain
% operation.
% 
% INPUTS:
% in	= input data in matlab structure format.  
% mode	= ppmrange('p') or Value('v').  In ppmrange mode, the DC offset is
%           calculated automatically over a defined ppm range in the 
%           spectrum.  In value mode, the user has to provide the value of 
%           the desired DC offset.
% var1	= If mode is 'p', then 'var1' is the range of ppm values 
%           [minppm maxppm] over which to calculate the DC offset.  In 
%           ppmrange mode, the 'var1' is optional, and takes a default 
%           value of [max(in.ppm)-0.5 max(in.ppm] AND 
%           [min(in.ppm) min(in.ppm)+0.5].  If mode is 'v', then 'var1' is 
%           the value of the dc offset correction that you wish to employ.
%
% OUTPUTS:
% out   = Output following DC correction.  

function [out,dcOffset]=op_dccorr(in,mode,var1);

if mode=='p'
    if nargin<3
        var1a=[max(in.ppm)-0.5 max(in.ppm)];
        var1b=[min(in.ppm) min(in.ppm)+0.5];
        dcOffset1=mean(in.specs(in.ppm<var1a(2) & in.ppm>var1a(1)));
        dcOffset2=mean(in.specs(in.ppm<var1b(2) & in.ppm>var1b(1)));
        dcOffset=mean([dcOffset1 dcOffset2]);
    else
        dcOffset=mean(in.specs(in.ppm<var1(2) & in.ppm>var1(1)));
    end
elseif mode=='v'
    dcOffset=var1;
end

specs=in.specs-dcOffset;

%convert back to time domain
%if the length of Fids is odd, then you have to do a circshift of one to
%make sure that you don't introduce a small frequency shift into the fids
%vector.
if mod(size(specs,in.dims.t),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fids=fft(fftshift(specs,in.dims.t),[],in.dims.t);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fids=fft(circshift(fftshift(specs,in.dims.t),1),[],in.dims.t);
end

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.freqranged=1;