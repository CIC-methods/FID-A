% sim_press_shaped.m
% Robin Simpson and Jamie Near, 2014.
% 
% USAGE:
% out = sim_press_shaped(n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dx,dy,Gx,Gy,phCyc1,phCyc2))
% 
% DESCRIPTION:
% This function simulates the PRESS experiment.  The excitation is
% simulated as an instantaneous rotation, and the refocusing pulse is
% simulated as a shaped rotation.
%
% This code enables the choice of the phase of the refocusing pulses.  This 
% enables phase cycling of the refocusing pulses by repeating simulations 
% with different editing pulse phases, which is necessary to remove phase 
% artefacts from the editing pulses.  A four step phase cycling scheme is typically
% sufficient, where both refocusing pulses are phase cycled by 0 and 90 degrees, and
% the phase are combined in the following way:
% 
% signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
% 
% where, in [X Y], X is the phase of the first refocusing pulse and Y is
% the phase of the second refocusing pulse
% 
% Finally, this code simulates the spectrum at a given point in space (x,y),
% given the values of the slice selection gradients (Gx, and Gy).  The pulse
% waveform is assumed to be the same for both refocusing pulses.  In order
% to fully simulate the MEGA-PRESS experiment, you have to run this
% simulation many times at various points in space (x,y), and then add
% together the resulting spectra.  
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% J         = matrix of coupling constants in [Hz] for spin system
% shifts    = vector of chemical shifts in [ppm] for spin system
% tau1      = echo time 1 in [ms].
% tau2      = echo time 2 in [ms].
% RF        = radiofrequency pulse array [N x 3].  Phase, Amplitude, Duration.
% tp        = RF pulse duration in [ms]
% peakB1    = peak B1 amplitude in [kHz]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = gradient strength for first selective refocusing pulse [G/cm]
% Gy        = gradient strength for second selective refocusing pulse [G/cm]
% phCycl    = initial phase of the first refocusing pulse in [degrees];
% phCycl2   = initial phase of the second refocusing pulse in [degrees];

function out = sim_press_shaped(n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dx,dy,Gx,Gy,phCyc1,phCyc2)
    
if tau1<tp/1000
    error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
end
if tau2<tp/1000
    error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
end

%Set water to centre
centreFreq=4.65;
sys.shifts=sys.shifts-centreFreq;

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN PULSE SEQUENCE************
d=sim_excite(H,'x');                                    %EXCITE
d=sim_evolve(d,H,tau1/2000);                            %Evolve by tau1/2
d=sim_shapedRF(d,H,RF,tp,180,90+phCyc1,Gx,dx);          %1st shaped 180 degree refocusing pulse
d=sim_evolve(d,H,(tau1+tau2)/2000);                     %Evolve by (tau1+tau2)/2
d=sim_shapedRF(d,H,RF,tp,180,90+phCyc2,Gy,dy);          %2nd shaped 180 degree refocusing pulse
d=sim_evolve(d,H,tau2/2000);                            %Evolve by tau2/2
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Fill in structure header fields:
out.seq='press';
out.te=tau1+tau2;
out.sim='shaped';

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








 
