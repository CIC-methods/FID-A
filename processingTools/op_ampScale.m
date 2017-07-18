% op_ampScale.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_ampScale(in,A);
% 
% DESCRIPTION:
% Scale the amplitude of a spectrum by factor A.
% 
% INPUTS:
% in    = input data in matlab structure format
% A     = Amplitude scaling factor.
%
% OUTPUTS:
% out   = Output following amplitude scaling.  

function out=op_ampScale(in,A);


out=in;
out.specs=in.specs*A;
out.fids=in.fids*A;
