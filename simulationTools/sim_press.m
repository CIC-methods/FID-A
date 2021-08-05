% sim_press.m
% Robin Simpson and Jamie Near, 2014.
%
% USAGE:
% out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2)
% 
% DESCRIPTION:
% This function simulates an ideal PRESS experiment with first echo time 
% "tau1" and a second echo time of "tau2".  The function calls the function
% 'sim_Hamiltonian.m' which produces the free evolution Hamiltonian and 
% rotation Hamiltonians for the specified spin system.
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% tau1      = Echo time in [ms] of first press Spin Echo
% tau2      = Echo time in [ms] of second press Spin Echo
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.

function out = sim_press(n,sw,Bfield,linewidth,sys,tau1,tau2)

%Set water to centre
centreFreq=4.65;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                            %EXCITE
d=sim_evolve(d,H,tau1/2000);                    %Evolve by tau1/2
d=sim_rotate(d,H,180,'y');                      %First 180 degree refocusing pulse about y' axis.
d=sim_evolve(d,H,(tau1+tau2)/2000);             %Evolve by (tau1+tau2)/2
d=sim_rotate(d,H,180,'y');                      %second 180 degree refocusing pulse about y' axis.
d=sim_evolve(d,H,tau2/2000);                    %Evolve by tau2/2
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);  %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='press';
out.te=(tau1+tau2);
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






 