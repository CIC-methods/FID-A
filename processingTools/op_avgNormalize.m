% op_avgNormalize.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_avgNormalize(in,scaleFactor);
% 
% DESCRIPTION:
% Divide by the number of averages.  This is typically performed after
% op_averaging.m to obtain a truely "averaged spectrum", since op_averaging.m
% adds the averages together but does not divide by the number of averages.
% 
% INPUTS:
% in            = input data in matlab structure format.
% scaleFactor	= (optional) Factor to divide by.  default = in.averages.

function out=op_avgNormalize(in,scaleFactor);

if in.flags.avgNormalized==1
    error('ERROR:  average normalization has already been performed!');
end

if nargin<2
    scaleFactor=in.averages;
end

out=in;
out.specs=in.specs/scaleFactor;
out.fids=in.fids/scaleFactor;
out.flags.avgNormalized=1;

