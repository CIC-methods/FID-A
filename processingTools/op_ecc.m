% op_ecc.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,outw]=op_ecc(in,inw);
% 
% DESCRIPTION:
% Perform an eddy current correction by estimating any non-linearity in the
% phase of the water unsuppressed data in the time domain and applying the
% appropriate correction to both the water suppressed and water unsuppressed 
% data.
% 
% INPUTS:
% in     = water suppressed input data in matlab structure format.
% inw    = water unsuppressed input data in matlab structure format.
%
% OUTPUTS:
% out    = Water suppressed output following eddy current correction  
% outw   = Water unsuppressed output following eddy current correction

function [out,outw]=op_ecc(in,inw);

if inw.dims.coils~=0 || inw.dims.averages~=0 || inw.dims.subSpecs~=0
    error('ERROR:  Must combine receivers, averages and subspecs prior to running ecc!! ABORTING!!');
end

%save the phase as a vector of hard numbers.
inph=double(phase(inw.fids));
plot(inw.t,inph);

%choose the part of the phase function that is most linear
tmin=input('input min t value: ');
tmax=input('input max t value: ');
figure;

%now fit a straight line to the linear part of the phase function
p=polyfit(inw.t(inw.t>tmin & inw.t<tmax),inph(inw.t>tmin & inw.t<tmax)',1)

%now fit a spline to approximate a smooth version of the phase function
pp=splinefit(inw.t,inph,150);

%now subtract the line from the spline to get the eddy current related
%phase offset:
ecphase=ppval(pp,inw.t)'-polyval(p,inw.t)';
sz=size(in.fids)
ecphase_rep=repmat(ecphase,[1 sz(2:end)]);
size(ecphase_rep)
figure;
plot(inw.t,ecphase);


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
figure;
plot(outw.t,phase(outw.fids));