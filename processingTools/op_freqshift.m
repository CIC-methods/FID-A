% op_freqshift.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_freqshift(in,f);
% 
% DESCRIPTION:
% Apply a frequency shift to the input spectrum by 'f' Hz.
% 
% INPUTS:
% in     = input data in matlab structure format.
% f      = frequency shift to apply (in Hz).
%
% OUTPUTS:
% out    = Output following frequency shift.  

function out=op_freqshift(in,f);


% if in.dims.coils>0
%     error('ERROR:  Can not operate on data with multilple coils!  ABORTING!!')
% end
% if in.dims.averages>0
%     error('ERROR:  Can not operate on data with multiple averages!  ABORTING!!');
% end
% if in.dims.subSpecs>0
%     error('ERROR:  Can not operate on data with multiple Subspecs!  ABORTING!!');
% end

t=repmat(in.t',[1 in.sz(2:end)]);

fids=in.fids.*exp(-1i*t*f*2*pi);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%plot(in1.ppm,combinedSpecs);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;