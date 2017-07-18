% op_movef0.m
% Jamie Near and Eric Krochmalnek, McGill University 2016.
% 
% USAGE:
% out=op_movef0(in,f);
% 
% DESCRIPTION:
% Change the centre frequency of the input spectrum by 'f' Hz.
% 
% INPUTS:
% in     = input data in matlab structure format.
% f      = frequency shift to apply (in Hz).
%
% OUTPUTS:
% out    = Output dataset following f0 shift.

function out=op_movef0(in,f);


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

%change the ppm scale as well
ppm=in.ppm+(f/(in.txfrq/1000000));

%plot(in1.ppm,combinedSpecs);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;
out.ppm=ppm;

%FILLING IN THE FLAGS
out.flags=in.flags;