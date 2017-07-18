% op_freqshiftSubspec.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_freqshiftSubspec(in,f);
% 
% DESCRIPTION:
% Apply a frequency shift to only one of the sub-spectra in a dataset.  This
% is used to minimize subtraction artefacts from MEGA_PRESS data, for
% instance.
% 
% INPUTS:
% in     = input data in matlab structure format.
% f      = frequency shift to apply to subspectrum (in Hz).
%
% OUTPUTS:
% out    = Output following frequency shifting of subspectrum.  

function out=op_freqshiftSubspec(in,f);

%this function is only meant to operate on scans with two subspecs.  It
%will apply a frequnency shift to the second subspectrum.  This is intended
%for use with edited spectroscopy sequences where small frequency drifts
%between edit-on and edit-off spectra can result in unwanted residual
%signals from uncoupled spins (Cr, Ch, etc.).

if in.dims.coils>0
    error('ERROR:  Can not operate on data with multilple coils!  ABORTING!!')
end
if in.dims.averages>0
    error('ERROR:  Can not operate on data with multiple averages!  ABORTING!!');
end
if in.dims.subSpecs==0
    error('ERROR:  Can not operate on data with no Subspecs!  ABORTING!!');
end
if in.sz(in.dims.subSpecs~=2)
    error('ERROR:  Input spectrum must have two subspecs!  ABORTING!!');
end

fids=in.fids;
fids(:,2)=fids(:,2).*exp(-1i*in.t'*f*2*pi);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%plot(in1.ppm,combinedSpecs);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;