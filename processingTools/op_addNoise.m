% op_addNoise.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,noise]=op_addNoise(in,sdnoise,noise);
% 
% DESCRIPTION:
% Add noise to a spectrum (useful for generating simulated data).  Normally
% distributed random noise is added to both the real and imaginary parts of
% the data.  Real and imaginary noise parts are uncorrelated.
% 
% INPUTS:
% in         = Input data in matlab structure format.
% sdnoise    = Standard deviation of the random noise to be added in the time domain.
% noise      = (optional)  Specific noise kernel to be added (if specified,
%               sdnoise variable is ignored).  
%
% OUTPUTS:
% out        = Output dataset with noise added.
% noise      = The noise vector that was added. 

function [out,noise]=op_addNoise(in,sdnoise,noise);

if nargin<3
    noise_real=sdnoise*randn(in.sz);
    noise_imag=sdnoise*randn(in.sz);
    noise=noise_real+(1i*noise_imag);
end

fids=in.fids + noise;

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;