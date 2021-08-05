%sim_megapress.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=sim_megapress(n,sw,Bfield,linewidth,sys,taus,refoc1Flip,refoc2Flip,editFlip)
% 
% DESCRIPTION:
% Simulate the MEGA-PRESS sequence with instantaneous localization and
% editing pulses.  Provides the ability to specify the flip angle of each
% refocusing pulse and editing pulse on each spin in the spin system.
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% taus      = pulse sequence timing vector:
%   taus(1)     = time in [ms] from 90 to 1st 180
%   taus(2)     = time in [ms] from 1st 180 to 1st edit pulse
%   taus(3)     = time in [ms] from 1st edit pulse to 2nd 180
%   taus(4)     = time in [ms] from 2nd 180 to 2nd edit pulse
%   taus(5)     = time in [ms] from 2nd edit pulse to ADC
% refoc1Flip= cell array of refoc1 flip angles for each spin in system
% refoc2Flip= cell array of refoc2 flip angles for each spin in system
% editFlip  = cell array of editing flip angles for each spin in system
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using MEGA-PRESS 
%             sequence.

function out = sim_megapress(n,sw,Bfield,linewidth,sys,taus,refoc1Flip,refoc2Flip,editFlip)

%Set water to centre
centreFreq=4.65;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                            %EXCITE
d=sim_evolve(d,H,taus(1)/1000);                      %Evolve by taus(1)
d=sim_rotate(d,H,refoc1Flip,'y');               %First refocusing pulse about y' axis.
d=sim_evolve(d,H,taus(2)/1000);                      %Evolve by taus(2)
d=sim_rotate(d,H,editFlip,'y');                 %First editing pulse about y' axis.
d=sim_evolve(d,H,taus(3)/1000);                      %Evolve by taus(3)
d=sim_rotate(d,H,refoc2Flip,'y');               %Second refocusing pulse about y' axis.
d=sim_evolve(d,H,taus(4)/1000);                      %Evolve by taus(4)
d=sim_rotate(d,H,editFlip,'y');                 %Second editing pulse about y' axis.
d=sim_evolve(d,H,taus(5)/1000);                      %Evolve by taus(5)
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);  %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='megapress';
out.te=sum(taus);
out.sim='ideal';

%Additional fields for compatibility with FID-A processing tools.
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










 