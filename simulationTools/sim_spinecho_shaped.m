% sim_spinecho_shaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% sim_spinecho_shaped(n,sw,Bfield,linewidth,sys,TE,RF,Tp,grad,pos,ph)
% 
% DESCRIPTION:
% This function simulates a localized spin-echo sequence with a shaped
% refocusing pulse.  It enables choice of the echo-time as well as the 
% choice of the phase refocusing pulse.  This allows phase cycling of the 
% refocusing pulses by repeating simulations with different pulse phases, 
% which is necessary to remove unwanted coherences from outside the volume 
% of interest.  For the refocusing pulse, a two step phase cycling scheme 
% is typically sufficient, where the refocusing pulse is phase cycled by 0 
% and 90 degrees the phase are combined by subtraction.
% 
% INPUTS:
% n          = number of points in fid/spectrum
% sw         = desired spectral width in [Hz]
% Bfield     = main magnetic field strength in [T]
% linewidth  = linewidth in [Hz]
% sys        = Metabolite spin system definition structure;
% TE         = Echo time in [ms]
% RF         = RF pulse definition structure for refoc pulse (obtain using 'io_loadRFwaveform.m');
% Tp         = duration of refocusing pulse in [ms]
% grad       = gradient strength for the selective refocusing pulse [G/cm]
% pos        = position offset in the direction corresponding to the refocusing pulse [cm]
% ph         = the phase of the refocusing pulse in [degrees];
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using spin-echo 
%             sequence.


function out = sim_spinecho_shaped(n,sw,Bfield,linewidth,sys,TE,RF,Tp,grad,pos,ph)

%Set water to centre
centreFreq=4.65;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delay by subtracting the pulse durations from TE;
delay=TE-Tp;
if delay<0
    error(['ERROR!  TE is too short.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                            %EXCITE
d=sim_evolve(d,H,delay/2000);                      %Evolve by delay/2
d=sim_shapedRF(d,H,RF,Tp,180,90+ph,pos,grad);   %shaped 180 degree refocusing pulse about y' axis.
d=sim_evolve(d,H,delay/2000);                      %Evolve by delay/2
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);  %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='spinecho';
out.te=TE;
out.sim='shaped';

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

