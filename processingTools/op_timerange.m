% op_timerange.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_timerange(in,tmin,tmax);
% 
% DESCRIPTION:
% Output only a specified frequency range of the input spectrum.
% 
% INPUTS:
% in         = input data in matlab structure format.
% tmin       = minimum extent of frequency range in ppm.
% tmax       = maximum extent of frequency range in ppm.
%
% OUTPUTS:
% out        = Output dataset following truncation in the time domain. 

function out=op_timerange(in,tmin,tmax);

%Take only the specified range of the FID
fids=in.fids(in.t>=tmin & in.t<tmax,:,:);

%Calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%calculate the size;
sz=size(fids);

%calculate the time scale
t=in.t(in.t>tmin & in.t<tmax);

%Now re-calculate t and ppm arrays using the calculated parameters:
f=[(-in.spectralwidth/2)+(in.spectralwidth/(2*sz(1))):...
    in.spectralwidth/(sz(1)):...
    (in.spectralwidth/2)-(in.spectralwidth/(2*sz(1)))];

ppm=-f/(in.Bo*42.577);
ppm=ppm+4.65;

%calculate the time scale
t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];



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
