% op_CSImakeK0fid.m
% Jamie Near & Sneha Senthil, Sunnybrook Research Insitute 2021.
% 
% USAGE:
% out=op_CSImakeK0fid(in,start,interval);
% 
% DESCRIPTION:
% Extract points with repeating periodic interval from a spatio-spectral 
% MRSI signal.  This is useful for Rosette MRSI data, when we want to extract 
% the FIDs from the centre of k-space.  
% 
% INPUTS:
% in        = input MRSI data in matlab structure format.
% start     = The index of the first point to extract.
% interval  = The interval between subsequently extracted points.
%
% OUTPUTS:
% out       = Extracted FID in single voxel FID-A data structure format.   

function out=op_CSImakeK0fid(in,start,interval);

%extract the desired points;
fids=in.data(start:interval:end,:,:,:);

%Get the new size:
sz=size(fids);

%remake the time vector:
t=in.adcTime(start:interval:end);

%Dwelltime and spectral width:
dwelltime = t(2)-t(1);
sw = 1/dwelltime;

%Calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%Now re-calculate ppm array using the calculated parameters:
f=[(-sw/2)+(sw/(2*sz(1))):...
    sw/(sz(1)):...
    (sw/2)-(sw/(2*sz(1)))];

ppm=-f/(in.Bo*42.577);
ppm=ppm+4.65;


%FILLING IN DATA STRUCTURE
out=in;
out=rmfield(out,'data');
out=rmfield(out,'adcTime');
out=rmfield(out,'adcDwellTime');
out=rmfield(out,'imageOrigin');
out=rmfield(out,'coordinates');
out=rmfield(out,'affineMatrix');
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.t=t;    
out.ppm=ppm;
out.dwelltime=dwelltime;
out.sw=1/dwelltime;


%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.isMRSI=0;
