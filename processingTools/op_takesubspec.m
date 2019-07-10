% op_takesubspec.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_takesubspec(in,index);
% 
% DESCRIPTION:
% Extract the subspectra with indices corresponding to the 'index' input
% array. 
% 
% INPUTS:
% in     = input data in matlab structure format.
% index  = vector indicating the indices of the subspectra you would like 
%          to extract.
%
% OUTPUTS:
% out    = Output dataset consisting of subspectra indices extracted from 
%          the input.

function out=op_takesubspec(in,index);

if in.flags.subtracted
    error('ERROR:  Subspectra have already been combined!  Aborting!');
end

if in.dims.subSpecs==0
    %Can't take subspec because there are none:
    error('ERROR:  There are no subspectra in this dataset!  Aborting!');
elseif in.dims.subSpecs==1
    %SHOULD NEVER HAPPEN (Time dimension should always be dim=1)
    error('ERROR:  dims.subSpecs==1.  This should never happen!  Aborting!');
elseif in.dims.subSpecs==2
    fids=in.fids(:,index);
elseif in.dims.subSpecs==3;
    fids=in.fids(:,:,index);
elseif in.dims.subSpecs==4;
    fids=in.fids(:,:,:,index);
elseif in.dims.subSpecs==5
    fids=in.fids(:,:,:,:,index);
end

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%change the dims variables
if in.dims.t>in.dims.subSpecs
    dims.t=in.dims.t-1;
else
    dims.t=in.dims.t;
end
if in.dims.coils>in.dims.subSpecs
    dims.coils=in.dims.coils-1;
else
    dims.coils=in.dims.coils;
end
if in.dims.averages>in.dims.subSpecs
    dims.averages=in.dims.averages-1;
else
    dims.averages=in.dims.averages;
end
if in.dims.extras>in.dims.subSpecs
    dims.extras=in.dims.extras-1;
else
    dims.extras=in.dims.extras;
end
dims.subSpecs=0;

%re-calculate the sz variable
sz=size(fids);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.dims=dims;
out.subspecs=1;
if out.dims.averages
    out.averages=out.sz(out.dims.averages);
else
    out.averages=1;
end
out.rawAverages=out.averages;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
