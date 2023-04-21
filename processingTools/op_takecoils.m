% op_takecoils.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_takecoils(in,index);
% 
% DESCRIPTION:
% Extract the rf coil channels with indices corresponding to the 'index' input
% array. 
% 
% INPUTS:
% in     = input data in matlab structure format.
% index  = vector indicating the indices of the rf channels you would like to extract.
%
% OUTPUTS:
% out    = Output dataset consisting of averages extracted from the input.

function out=op_takecoils(in,index);

if in.dims.coils==0
    %Can't take average because there are none:
    error('ERROR:  There aren''t multiple coils in this dataset!  Aborting!');
elseif in.dims.coils==1
    %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
    error('ERROR:  dims.coils==1.  This should never happen!  Aborting!');
elseif in.dims.coils==2
    fids=in.fids(:,index,:,:,:);
elseif in.dims.coils==3;
    fids=in.fids(:,:,index,:,:);
elseif in.dims.coils==4;
    fids=in.fids(:,:,:,index,:);
elseif in.dims.coils==5
    fids=in.fids(:,:,:,:,index);
end

fids=squeeze(fids);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
dims.t=in.dims.t;
if length(index)==1
    dims.coils=0;
    if in.dims.averages > in.dims.coils
        dims.averages=in.dims.averages-1;
    else
        dims.averages=in.dims.averages;
    end
    if in.dims.subSpecs > in.dims.coils
        dims.subSpecs=in.dims.subSpecs-1;
    else
        dims.subSpecs=in.dims.subSpecs;
    end
    if in.dims.extras > in.dims.coils
        dims.extras=in.dims.extras-1
    else
        dims.extras=in.dims.extras;
    end
elseif length(index)>1
    dims.coils=in.dims.coils;
    dims.averages=in.dims.averages;
    dims.subSpecs=in.dims.subSpecs;
    dims.extras=in.dims.extras;
end

%re-calculate the sz variable
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
if length(index)<2
    out.flags.addedrcvrs=1;
end

