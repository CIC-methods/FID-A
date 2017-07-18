% op_addphaseSubspec.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_addphaseSubspec(in,ph);
% 
% DESCRIPTION:
% Add zero order phase to one of the subspectra in a dataset.  For example,
% the edit-on spectrum of a mega-press acquisition.
% 
% INPUTS:
% in    = Input spectrum in matlab structure format.
% ph    = Phase (in degrees) to add to the second subspectrum.
%
% OUTPUTS:
% out   = Output dataset with phase adjusted subspectrum. 


function out=op_addphaseSubspec(in,ph);

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
fids(:,2)=fids(:,2).*exp(1i*ph*pi/180);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%plot(in1.ppm,combinedSpecs);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;