% op_takeextras.m
% Jamie Near, McGill University 2018.
% 
% USAGE:
% out=op_takeextras(in,index);
% 
% DESCRIPTION:
% Extract the extras with indices corresponding to the 'index' input
% array. 
% 
% INPUTS:
% in     = input data in matlab structure format.
% index  = vector indicating the indices of the extras you would like to extract.
%
% OUTPUTS:
% out    = Output dataset consisting of extras extracted from the input.

function out=op_takeextras(in,index);

if in.dims.extras==0
    %Can't take extras because there are none:
    error('ERROR:  There are no averages in this dataset!  Aborting!');
elseif in.dims.extras==1
    %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
    error('ERROR:  dims.subSpecs==1.  This should never happen!  Aborting!');
elseif in.dims.extras==2
    fids=in.fids(:,index,:,:,:);
elseif in.dims.extras==3;
    fids=in.fids(:,:,index,:,:);
elseif in.dims.extras==4;
    fids=in.fids(:,:,:,index,:);
elseif in.dims.extras==5
    fids=in.fids(:,:,:,:,index);
end

fids=squeeze(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
dims.t=in.dims.t;
if length(index)==1
    dims.extras=0;
    if in.dims.coils > in.dims.extras
        dims.coils=in.dims.coils-1;
    else
        dims.coils=in.dims.coils;
    end
    if in.dims.subSpecs > in.dims.extras
        dims.subSpecs=in.dims.subSpecs-1;
    else
        dims.subSpecs=in.dims.subSpecs;
    end
    if in.dims.averages > in.dims.extras
        dims.averages=in.dims.averages-1;
    else
        dims.averages=in.dims.averages;
    end
elseif length(index)>1
    dims.coils=in.dims.coils;
    dims.averages=in.dims.averages;
    dims.subSpecs=in.dims.subSpecs;
    dims.extras=in.dims.extras;
end

%re-calculate the sz variable
sz=size(fids);

%change the number of extras
if dims.extras==0
    extras=1;
else
    extras=sz(dims.extras);
end

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
