% op_fddccorr.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_fddccorr(in,npts);
% 
% DESCRIPTION:
% Correct and DC offset in the frequency domain.  This is equivalent to a
% vertical shift in the frequency domain.  The required vertical shift is
% calculated by taking the average of the first and last "NPTS" points in
% the frequency domain.  This requires that those points are in the noise
% floor.  
% 
% INPUTS:
% in     = input data in matlab structure format.
% npts   = number of points at both edges of the freqeuncy domain that will
%         be used for estimation of the DC offset of the spectrum.
%
% OUTPUTS:
% out    = Output following time-domain DC offset correction.  

function out=op_fddccorr(in,npts);


%Find the frequency domain vertical offset to correct.
tails=([in.specs(1:npts); in.specs(end-npts+1:end)]);
offset=mean(tails);

%subtract the offset (real and imaginary) from specs.
specs=in.specs;
specs=specs-offset;

%recalculate fids using fft
fids=fft(ifftshift(specs));


    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;

