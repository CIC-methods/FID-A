% op_ecc_klose.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,outw]=op_ecc_klose(in,inw);
% 
% DESCRIPTION:
% Perform an eddy current correction by subtracting the phase of the water
% FID, as suggested by Klose et al. 
% 
% INPUTS:
% in     = water suppressed input data in matlab structure format.
% inw    = water unsuppressed input data in matlab structure format.
%
% OUTPUTS:
% out    = Water suppressed output following eddy current correction  
% outw   = Water unsuppressed output following eddy current correction

function [out,outw]=op_ecc_klose(in,inw)

if inw.dims.coils~=0 || inw.dims.averages~=0 || inw.dims.subSpecs~=0
    error('ERROR:  Must combine receivers, averages and subspecs prior to running ecc!! ABORTING!!');
end

%save the phase as a vector of hard numbers.
inph=double(phase(inw.fids));
% plot(inw.t,inph);

%now subtract the line from the spline to get the eddy current related
%phase offset:
ecphase=inph;
sz=size(in.fids);
ecphase_rep=repmat(ecphase,[1 sz(2:end)]);
size(ecphase_rep);


%Now apply the eddy current correction to both the water suppressed and the
%water unsuppressed data:
out=in;
out.fids=out.fids.*exp(1i*-ecphase_rep);
out.specs=fftshift(ifft(out.fids,[],1),1);
size(ecphase);
ecphase(1);
out=op_addphase(out,180*ecphase_rep(1)/pi);

outw=inw;
outw.fids=outw.fids.*exp(1i*-ecphase);
outw.specs=fftshift(ifft(outw.fids,[],1),1);
outw=op_addphase(outw,180*ecphase(1)/pi);
% figure;
% plot(outw.t,phase(outw.fids));