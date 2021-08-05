%sim_readout.m
%Jamie Near, 2014.
%Speed-up modifications by Chathura Kumaragamage, Feb 2018.
%
% USAGE:
% [out,d_out] = sim_readout(d_in,H,n,sw,linewidth,rcvPhase,shape)
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
%                'LG' = Lorentz-Gauss (50% mixture)
%
% OUTPUTS:
% out       = simulated spectrum resulting from readout.
% d_out     = output density matrix following readout.

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
Bfield=H(1).Bfield;
if strcmp(shape,'L') || strcmp(shape,'l')
    t2 = 1/(pi*linewidth);
elseif strcmp(shape,'G') || strcmp(shape,'g')
    thalf=log(0.5)/(pi*0.5*linewidth);
    sigma=sqrt((thalf^2)/(-2*log(0.5)));
elseif strcmp(shape,'LG') || strcmp(shape,'lg')
    t2 = 1/(pi*linewidth);
    thalf=log(0.5)/(pi*0.5*linewidth);
    sigma=sqrt((thalf^2)/(-2*log(0.5)));
    R=0.5;
else 
    error('ERROR:  Shape not recognized!  ABORTING!');
end

%READOUT CODE COURTESY OF SAAD JBABDI, FMRIB 2011!!
%FURTHER SPEEDUP IMPLEMENTED BY CHATHURA KUMARAGAMAGE, MCGILL 2018!!
%LOOPING THROUGH ARRAYED SPIN SYSTEMS.  JAMIE NEAR, MCGILL 2018!!
l=0:points-1;
if strcmp(shape,'L') || strcmp(shape,'l')
    decay=exp(-((l)*deltat)/t2);
elseif strcmp(shape,'G') || strcmp(shape,'g')
    decay = exp(-(((l)*deltat).^2)/(2*(sigma^2)));
elseif strcmp(shape,'LG') || strcmp(shape,'lg')
    decay = (R*exp(-((l)*deltat)/t2))+((1-R)*exp(-(((l)*deltat)^2)/(2*(sigma^2))));
end

for j=1:length(H) %JN - Loop through the different parts of the spin-system:
    M=H(j).HAB;
    [U,D]=eig(M);D=diag(D);
    val=2^(2-H(j).nspins);
    k=1;

    Fxy =H(j).Fx+1i*H(j).Fy;
    phase_comp=exp(1i*rcvPhase*pi/180);
    
    while k<(points+1);
        d1=diag(exp(-1i*(k-1)*D*deltat));
        d=U*d1*U';
        out_parts{j}.fids(k)=trace(d*d_in{j}*d'*(Fxy)*phase_comp);
        k = k+1;
    end

    out_parts{j}.fids = val*(out_parts{j}.fids.*decay);
        
end

d_out=sim_evolve(d_in,H,points*deltat);

%Combine the output FIDs from the different parts of the spin-system:
out.fids=zeros(size(out_parts{length(H)}.fids));
for m=1:length(H)
    out.fids=out.fids+out_parts{m}.fids;
end

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
out.txfrq=out.Bo*42577000;  %assumes proton.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
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
out.flags.isFourSteps=0;

