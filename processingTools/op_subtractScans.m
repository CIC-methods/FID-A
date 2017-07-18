% op_subtractScans.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_subtractScans(in1,in2);
% 
% DESCRIPTION:
% Subtract input 2 from input 1;
% 
% INPUTS:
% in1        = 1st input data in matlab structure format.
% in2        = 2nd input data in matlab structure format.
%
% OUTPUTS:
% out        = Output dataset following subtraction of in2 from in1.

function out=op_subtractScans(in1,in2);

fids=in1.fids-in2.fids;

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in1.dims.t),in1.dims.t);

%FILLING IN DATA STRUCTURES
out=in1;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in1.flags;