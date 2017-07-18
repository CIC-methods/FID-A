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
sz=size(fids);

fids=conj(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%Calculate t and ppm arrays using the calculated parameters:
f=[(-in.spectralwidth/2)+(in.spectralwidth/(2*sz(1))):in.spectralwidth/(sz(1)):(in.spectralwidth/2)-(in.spectralwidth/(2*sz(1)))];
%ppm=-f/(in.Bo*42.577);
ppm=-f/(3*42.577);
ppm=ppm+4.65;

%t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];
t=[in.dwelltime:in.dwelltime:sz(1)*in.dwelltime];

    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
