%sim_steam.m
%Jamie Near, McGill University 2014; Dana Goerzen, McGill University 2021
%
% USAGE:
% out = sim_steam(n,sw,Bfield,linewidth,sys,te,tm)
% 
% DESCRIPTION:
% Simulate the STEAM sequence using ideal (instantaneous) RF pulses.  To 
% remove unwanted coherences, coherence order filtering is employed 
% THIS CODE IS NOT TESTED.  RESULTS MAY NOT BE ACCURATE!!
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% te        = echo time in [ms]
% tm        = mixing time in [ms]
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using steam
%             sequence.

function out = sim_steam(n,sw,Bfield,linewidth,sys,te,tm)

%Set water to centre
centreFreq=4.65;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-4.65;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

out=struct();
    %BEGIN PULSE SEQUENCE************
    d=sim_excite(d,H,'x');                             %EXCITE
    d=sim_COF(H,d,1);                                  %Select coherence order 1 
    d=sim_evolve(d,H,te/2000);                         %Evolve by te/2
    d=sim_rotate(d,H,-90,'x');                         %Second 90 degree pulse about x' axis.
    d=sim_COF(H,d,0);                                  %Select coherence order 0
    d=sim_evolve(d,H,tm/1000);                         %Evolve by TM delay
    d=sim_rotate(d,H,90,'x');                          %Second 90 degree pulse about x' axis.
    d=sim_COF(H,d,-1);                                 %Select coherence order -1
    d=sim_evolve(d,H,te/2000);                         %Evolve by te/2
    [out,dout]=sim_readout(d,H,n,sw,linewidth,90);     %Readout along +y' (90 degree phase);
%END PULSE SEQUENCE**************
%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='steam';
out.te=te;
out.tm=tm;
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
end
