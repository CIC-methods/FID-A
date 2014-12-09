%op_takesubspec.m
%Jamie Near, McGill University 2014.
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
% index  = vector indicating the indices of the subspectra you would like to extract.

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
dims.t=in.dims.t;
dims.coils=in.dims.coils;
dims.averages=in.dims.averages;
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

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.subtracted=1;