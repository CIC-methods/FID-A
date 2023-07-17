% op_complexConj.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_complexConj(in)
% 
% DESCRIPTION:
% take the complex conjugate of the data;
% 
% INPUTS:
% in	= Input data in matlab structure format.
%
% OUTPUTS:
% out   = Output following conjugation.  

function out=op_complexConj(in);

fids=in.fids;
fids=conj(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];
t=[in.dwelltime:in.dwelltime:sz(1)*in.dwelltime];

    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
