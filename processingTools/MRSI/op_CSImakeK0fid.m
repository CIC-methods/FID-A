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
% out       = Extracted FID.   

function out=op_CSImakeK0fid(in,start,interval);

%extract the desired points;
fids=in.data(start:interval:end,:,:,:);

%recalculate the size;
sz=size(fids);

%remake the time vector:
t=in.adcTime(start:interval:end);

%dwelltime
dwelltime=t(2)-t(1);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.sz=sz;
out.t=t;    
out.dwelltime=dwelltime;
out.sw=1/dwelltime;


%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
