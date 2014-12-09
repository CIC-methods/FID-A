%sim_readout.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_readout(d_in,H,n,sw,linewidth,rcvPhase,shape)
% 
% DESCRIPTION:
% This function simulates an ADC readout of the transverse magnetization
% during the free evolution of the spin system under the effects of
% chemical shift and scalar coupling.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% n         = number of readout points
% sw        = spectral width [Hz]
% linewidth = full width at half maximum of spectral peaks [Hz] 
% rcvPhase  = receiver phase [degrees].  Optional.  Default = 0 (corresponds to x'-axis readout);
% shape     = line broadening function.  Optional,
%                'L' = lorentzian (default) 
%                'G' = gaussian 

function [out,d_out] = sim_readout(d_in,H,n,sw,linewidth,rcvPhase,shape);

if nargin<7
    shape='L'; %Lorentzian lineshape by default;
    if nargin<6
        rcvPhase=0;
    end
end

%initialize parameters;
deltat = 1/sw;
points = n;
out.fids = zeros(1,points);
Bfield=H.Bfield;
if shape=='L' || shape =='l'
    t2 = 1/(pi*linewidth);
elseif shape=='G' || shape=='g'
    thalf=log(0.5)/(pi*0.5*linewidth);
    sigma=sqrt((thalf^2)/(-2*log(0.5)));
else 
    error('ERROR:  Shape not recognized!  ABORTING!');
end

%READOUT CODE COURTESY OF SAAD JBABDI, FMRIB 2011!!
M=H.HAB;
[U,D]=eig(M);D=diag(D);
val=2^(2-H.nspins);
k=1;
while k<(points+1);
    d1=diag(exp(-1i*(k-1)*D*deltat));
    d2=diag(exp(1i*(k-1)*D*deltat));
    if shape=='L' || shape=='l'
        out.fids(k)=val*exp(-((k-1)*deltat)/t2)*trace(U*d1*U'*d_in*U*d2*U'*(H.Fx+1i*H.Fy)*exp(1i*rcvPhase*pi/180));
    elseif shape=='G' || shape=='g'
        out.fids(k)=val*exp(-(((k-1)*deltat)^2)/(2*(sigma^2)))*trace(U*d1*U'*d_in*U*d2*U'*(H.Fx+1i*H.Fy)*exp(1i*rcvPhase*pi/180));
    end  
    k = k+1;
end
d_out=sim_evolve(d_in,H,points*deltat);
out.t=[0:deltat:deltat*(points-1)];
freq=[(-sw/2)+(sw/(2*points)):sw/(points):(sw/2)-(sw/(2*points))];
ppm=-freq/(Bfield*42.577);
ppm=ppm+4.65;
out.ppm=ppm;

out.fids=out.fids';
out.specs=fftshift(ifft(out.fids));


%Add a few additional bits to ouput structure so that it matches the
%configuration of experimental data structures (loaded with jn_loadspec.m)
out.spectralwidth=sw;
out.dwelltime=deltat;
out.n=n;
out.linewidth=linewidth;
out.Bo=Bfield;
out.txfrq=out.Bo*42.577;
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.avgNormalized=1;
out.flags.isISIS=0;

