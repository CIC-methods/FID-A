% op_dccorr.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_dccorr(in,mode,var1);
% 
% DESCRIPTION:
% Do a DC Correction on the data.  This method is a frequency domain
% operation.
% 
% INPUTS:
% in	= input data in matlab structure format.  
% mode	= Point('p') or Value('v').  In point mode, the DC offset is
%           calculated automatically at a specific point in the spectrum.  In value
%           mode, the user has to provide the value of the desired DC offset.
% var1	= If mode is 'p', then 'var' is the index of the spectral point that you
%           wish to use to calculate the DC offset.  If mode is 'v', then
%           'var' is the value of the dc offset correction that you wish to
%           employ.
%
% OUTPUTS:
% out   = Output following DC correction.  

function out=op_dccorr(in,mode,var1);

if mode=='p'
    infilt=op_filter(in,5);
    dcOffset=infilt.specs(var1);
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