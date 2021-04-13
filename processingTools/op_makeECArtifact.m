% op_makeECArtifact.m
% Jamie Near, McGill University 2019.
% 
% USAGE:
% [out,fDrift]=op_makeECArtifact(in,A,tc);
% 
% DESCRIPTION:
% Add some synthetic eddy current artefacts to an MRS dataset.  The eddy
% current artefact is a damped exponential of the form EC = A*exp(-t/tc), 
% by which the frequency of the FID is modulated.  A is the amplitude 
% (in Hz) and tc is the time constant (in seconds).  
% 
% INPUTS:
% in         = input data in matlab structure format.
% A          = Amplitude (in Hz) of the eddy current artefact in the
%              time domain
% tc         = Time constant (in seconds) of the exponentially decaying
%              phase artefact in the time domain.  
%
% OUTPUTS:
% out        = Output dataset with eddy current artefact added.
% fFunc      = The frequency modulation that was added to the FID.

function [out,fFunc]=op_makeECArtifact(in,A,tc);
%
% if in.dims.coils>0
%     error('ERROR:  Can not operate on data with multilple coils!  ABORTING!!')
% end
% if in.dims.averages>0
%     error('ERROR:  Can not operate on data with multiple averages!  ABORTING!!');
% end
% if in.dims.subSpecs>0
%     error('ERROR:  Can not operate on data with multiple Subspecs!  ABORTING!!');
% end

%Make the eddy-current induced frequency modulation function:
fFunc=A*exp(-in.t/tc);

%Replicate the time and frequency modulation functions to the correct dimensions
%so that we can work with them:
tRep=repmat(in.t',[1 in.sz(2:end)]);
fRep=repmat(fFunc',[1 in.sz(2:end)]);

%apply the eddy current artefact:
fids=in.fids.*exp(-1i*tRep.*fRep*2*pi);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%plot(in1.ppm,combinedSpecs);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;




