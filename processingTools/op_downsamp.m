% op_downsamp.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_downsamp(in,dsFactor);
% 
% DESCRIPTION:  Change the time domain sampling rate of a spectrum by a
% factor of 'dsFactor'.  Nearest neighbour interpolation is performed by
% default.
% 
% INPUTS:
% in         = input data in matlab structure format.
% dsFactor   = factor by which to divide the sampling rate of the fid.
%
% OUTPUTS:
% out        = Output following downsampling 

function out=op_downsamp(in,dsFactor);

if length(in.sz)>2
    error('ERROR:  must combine averages, subspecs and coils first!');
end


%Add zeros using MATLAB array zeropadding function;
fids=resample(double(in.fids),1,dsFactor,0);

%Calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%recalculate the sz vector
sz=size(fids);

%recalculate the dwell time and spectral width
dwelltime=in.dwelltime*dsFactor;
spectralwidth=in.spectralwidth/dsFactor;


%Now re-calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):...
    spectralwidth/(sz(1)):...
    (spectralwidth/2)-(spectralwidth/(2*sz(1)))];

ppm=-f/(in.Bo*42.577);
ppm=ppm+4.65;

t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;
out.dwelltime=dwelltime;
out.spectralwidth=spectralwidth;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.downsampled=1;
