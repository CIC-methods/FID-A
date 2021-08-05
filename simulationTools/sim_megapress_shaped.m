% sim_megapress_shaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% sim_megapress_shaped(n,sw,Bfield,linewidth,taus,sys,editPulse,editTp,editPh1,editPh2,refPulse,refTp,Gx,Gy,dx,dy,refPh1,refPh2)
%
% DESCRIPTION:
% This function simulates the MEGA-PRESS sequence with shaped
% localization pulses and shaped editing pulses.  Enables choice of the 
% timings of all of the rf pulses as well as the choice of the
% phase of both the editing pulse and the refocusing pulses.  This allows 
% phase cycling of the editing and refocusing pulses by repeating 
% simulations with different editing pulse phases, which is necessary to remove phase 
% artefacts from the editing pulses.  For the editing pulses, an eight step 
% phase cycling scheme is typically sufficient, where by the first editing pulse is cycled by 
% 0 and 90 degrees, and the second editing pulse is cycled by 0,90,180, and 270 degrees, and all
% phase cycles should be added together to remove unwanted coherences.  For
% the refocusing pulses, a four step phase cycling scheme is typically
% sufficient, where both refocusing pulses are phase cycled by 0 and 90 degrees, and
% the phase are combined in the following way:
% 
% signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
% 
% where, in [X Y], X is the phase of the first refocusing pulse and Y is
% the phase of the second refocusing pulse
% 
% Note that this code only simulates one subspectrum at a time (edit-on or
% edit-off).  The difference spectrum can be obtained by simulating one of
% each, and then subtracting.
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% taus(1)     = time in [ms] from 90 to 1st 180
% taus(2)     = time in [ms] from 1st 180 to 1st edit pulse
% taus(3)     = time in [ms] from 1st edit pulse to 2nd 180
% taus(4)     = time in [ms] from 2nd 180 to 2nd edit pulse
% taus(5)     = time in [ms] from 2nd edit pulse to ADC
%               FOR MEGA-PRESS on SIEMENS SYSTEM:  taus=[4.545,12.7025,21.7975,12.7025,17.2526];
% sys        = Metabolite spin system definition structure;
% editPulse  = RF pulse definition structure for editing pulses (obtain using 'io_loadRFwaveform.m')
% editTp     = duration of editing pulse in [ms];
% editPh1    = the phase of the first editing pulse in [degrees];
% editPh2    = the phase of the second editing pulse in [degrees];
% refPulse   = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% refTp      = duration of refocusing pulse in [ms]
% Gx         = gradient strength for first selective refocusing pulse [G/cm]
% Gy         = gradient strength for second selective refocusing pulse [G/cm]
% dx         = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy         = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% refPh1     = the phase of the first refocusing pulse in [degrees];
% refPh2     = the phase of the second refocusing pulse in [degrees];
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using MEGA-PRESS 
%             sequence.

function out = sim_megapress_shaped(n,sw,Bfield,linewidth,taus,sys,editPulse,editTp,editPh1,editPh2,refPulse,refTp,Gx,Gy,dx,dy,refPh1,refPh2)

%Set 3ppm GABA resonance to centre
centreFreq=3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse durations from the taus
%vector;
delays=zeros(size(taus));
delays(1)=taus(1)-(refTp/2);
delays(2)=taus(2)-((refTp+editTp)/2);
delays(3)=taus(3)-((editTp+refTp)/2);
delays(4)=taus(4)-((refTp+editTp)/2);
delays(5)=taus(5)-(editTp/2);
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                     %EXCITE
d=sim_evolve(d,H,delays(1)/1000);                          %Evolve by delays(1)
d=sim_shapedRF(d,H,refPulse,refTp,180,90+refPh1,dx,Gx);  %1st shaped 180 degree refocusing pulse
d=sim_evolve(d,H,delays(2)/1000);                          %Evolve by delays(2)
d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh1);     %1st shaped editing pulse rotation
d=sim_evolve(d,H,delays(3)/1000);                          %Evolve by delays(3)
d=sim_shapedRF(d,H,refPulse,refTp,180,90+refPh2,dy,Gy);  %2nd shaped 180 degree refocusing pulse
d=sim_evolve(d,H,delays(4)/1000);                          %Evolve by delays(4)
d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh2);     %2nd shaped editing pulse rotation
d=sim_evolve(d,H,delays(5)/1000);                          %Evolve by delays(5)
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);           %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='megapress';
out.te=sum(taus);
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


