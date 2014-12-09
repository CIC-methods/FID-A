%op_takeaverages.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=op_takeaverages(in,index);
% 
% DESCRIPTION:
% Extract the averages with indices corresponding to the 'index' input
% array. 
% 
% INPUTS:
% in     = input data in matlab structure format.
% index  = vector indicating the indices of the averages you would like to extract.

function out=op_takeaverages(in,index);

if in.dims.averages==0
    %Can't take average because there are none:
    error('ERROR:  There are no averages in this dataset!  Aborting!');
elseif in.dims.averages==1
    %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
    error('ERROR:  dims.subSpecs==1.  This should never happen!  Aborting!');
elseif in.dims.averages==2
    fids=in.fids(:,index,:,:,:);
elseif in.dims.averages==3;
    fids=in.fids(:,:,index,:,:);
elseif in.dims.averages==4;
    fids=in.fids(:,:,:,index,:);
elseif in.dims.averages==5
    fids=in.fids(:,:,:,:,index);
end

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
dims.t=in.dims.t;
dims.coils=in.dims.coils;
if length(index)==1
    dims.averages=0;
elseif length(index)>1
    dims.averages=in.dims.averages;
end
dims.subSpecs=in.dims.subSpecs;

%re-calculate the sz variable
sz=size(fids);

%change the number of averages
if dims.averages==0
    averages=1;
else
    averages=sz(dims.averages);
end

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;
out.averages=averages;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
