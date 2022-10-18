% sim_steam_shaped.m
% Jamie Near, Sunnybrook Research Institute 2022.
% 
% USAGE:
% out = sim_steam_shaped(n,sw,Bfield,linewidth,sys,TE,TM,RF,tp,dx,dy,Gx,Gy,flipAngle,centreFreq)
% 
% DESCRIPTION:
% This function simulates the STEAM experiment.  The initial excitation is
% simulated as an instantaneous rotation, and the subsequent 90 degree 
% pulses are simulated as shaped rotations.
%
% It employs coherence selection to only include desired coherence orders
% of the spin system being simulated. For the STEAM experiment, the desired
% coherence after the excitation pulse is +1, 0 after the first refocusing
% pulse and -1 after the final refocusing pulse.
% 
% Finally, this code simulates the spectrum at a given point in space (x,y),
% given the values of the slice selection gradients (Gx, and Gy).  The pulse
% waveform is assumed to be the same for both refocusing pulses.  In order
% to fully simulate the STEAM experiment, you have to run this
% simulation many times at various points in space (x,y), and then add
% together the resulting spectra.  
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% TE        = echo time in [ms].
% TM        = mixing time in [ms].
% RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = gradient strength for first selective refocusing pulse [G/cm]
% Gy        = gradient strength for second selective refocusing pulse [G/cm]
% flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 90 deg)
% centreFreq= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.

function out = sim_steam_shaped(n,sw,Bfield,linewidth,sys,TE,TM,RF,tp,dx,dy,Gx,Gy,flipAngle,centreFreq)

if nargin<15
    centreFreq=2.3;
    if nargin<14
        flipAngle=90;
    end
end
   
%In the steam sequence, it can be common to use an asymmetric RF pulse for the
%90-degree pulse waveform.  In this case, it is conventional for the 2nd
%and 3rd rf pulses to be time-reversed versions of eachother, with the 2nd
%pulse (RF1) being a max-phase pulse, and the 3rd pulse (RF2) being a 
%min-phase pulse (i.e. the long tails of both pulse occurring during the TM 
%period so that the TE is minimized).  Here, use the asymmetry factor of 
%the pulse to make sure that the 2nd and 3rd pulses are max-phase and 
%min-phase, respectively:
if RF.rfCentre>0.5
    RF1=rf_timeReverse(RF);
    RF2=RF;
else
    RF1=RF;
    RF2=rf_timeReverse(RF);
end

%Check that the TE and TM values are not too short
if TE<(RF1.rfCentre*tp*2)
    error('ERROR:  TE cannot be less than duration of RF pulse! ABORTING!!');
end
if TM<(RF2.rfCentre*tp*2)
    error('ERROR:  Echo-time 2 cannot be less than duration of RF pulse! ABORTING!!');
end

%Set centre frequency
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse duration from TE and TM;
delays=zeros(2);
delays(1)=TE-(RF1.rfCentre*tp*2);
delays(2)=TM-(RF2.rfCentre*tp*2);
if sum(delays<0)
    error(['ERROR! The following timings are too short: ' num2str(find(delays<0)) '.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                    %EXCITE
d=sim_COF(H,d,1);                                       %Keep only +1-order coherences
d=sim_evolve(d,H,delays(1)/2000);                         %Evolve by delays(1)/2
d=sim_gradSpoil(d,H,[Gx,0,0],[dx,dy,0],tp*RF1.rfCentre);     %Prewind gradient for 2nd 90 degree pulse (Not sure why, but this only works when Gx is positive.  Intiutively, Gx amplitude should be the opposite of the slice select gradient (i.e. -Gx), but this does not seem to work).  
d=sim_shapedRF(d,H,RF1,tp,flipAngle,90,dx,Gx);             %1st shaped 90 degree selective pulse
d=sim_COF(H,d,0);                                       %Keep only 0-order coherences
d=sim_evolve(d,H,(delays(2))/1000);                       %Evolve by delays(2)
d=sim_shapedRF(d,H,RF2,tp,flipAngle,90,dy,Gy);           %2nd shaped 90 degree selective pulse
d=sim_gradSpoil(d,H,[0,Gy,0],[dx,dy,0],tp*RF1.rfCentre);    %Rewind gradient for 3rd 90 degree pulse (Not sure why, but this only works when Gy is positive.  Intiutively, Gx amplitude should be the opposite of the slice select gradient (i.e. -Gy), but this does not seem to work).
d=sim_COF(H,d,-1);                                      %Keep only -1 coherences
d=sim_evolve(d,H,delays(1)/2000);                         %Evolve by delays(1)/2
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);          %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='steam';
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


