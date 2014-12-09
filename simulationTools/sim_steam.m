%sim_steam.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out = sim_steam(n,sw,Bfield,linewidth,sys,te,tm,spoilAngle)
% 
% DESCRIPTION:
% Simulate the STEAM sequence using a spoiler gradient to get rid of
% unwanted coherences.  Simulation should be run multiple times with
% different spoil angles to get rid of all unwanted coherences.  THIS CODE
% HAS NOT BEEN TESTED.  NOT SURE IF IT PRODUCES THE CORRECT SPECTRA.
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% te        = echo time in [s]
% tm        = mixing time in [s]

function out = sim_steam(n,sw,Bfield,linewidth,sys,te,tm,spoilAngle)

%Set water to centre
sys.shifts=sys.shifts-4.65;

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN PULSE SEQUENCE************
d=sim_excite(H,'x');                            %EXCITE
d=sim_evolve(d,H,te/2);                         %Evolve by te/2
d=sim_spoil(d,H,spoilAngle);                    %Apply first te/2 spoiling;
d=sim_rotate(d,H,90,'x');                       %Second 90 degree pulse about x' axis.
d=sim_evolve(d,H,tm);                           %Evolve by TM delay
d=sim_spoil(d,H,spoilAngle);                    %Apply TM spoiling
d=sim_rotate(d,H,90,'x');                       %Final 90 degree pulse about x' axis.
d=sim_evolve(d,H,te/2);                         %Evolve by te/2
d=sim_spoil(d,H,spoilAngle);                    %Apply final te/2 spoiling
[out,dout]=sim_readout(d,H,n,sw,linewidth,270); %Readout along -y' (270 degree phase);
%END PULSE SEQUENCE**************

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














 