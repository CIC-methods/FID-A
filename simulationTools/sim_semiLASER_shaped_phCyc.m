% 2021 EDIT: This function is now deprecated. The current version of the 
% shaped PRESS simulation is sim_press_shaped_coFilt.m. The current version employs
% coherence selection rather than phase cycling to null undesired signal.
    
% sim_semiLASER_shaped_phCyc.m
% Dana Goerzen (McGill University, 2019).
% 
% USAGE:
%  out = sim_semiLASER_shaped_phCyc(n,sw,Bfield,linewidth,sys,te,RF,tp,dx,dy,Gx,Gy,ph1,ph2,ph3,ph4,flipAngle,centreFreq)
% 
% DESCRIPTION:
% This function simulates the semi-LASER experiment as described by Oz et al. (2018).
% The excitation is simulated as an instantaneous rotation, and the two pairs of slice
% selective refocusing adiabatic pulse are simulated as shaped RF pulses.

% This code enables the choice of the phase of the refocusing pulses.
% A four step phase cycling scheme is typically sufficient, where both 
% pairs of refocusing pulses are phase cycled by  0 and 90 degrees, 
% and the phases are combined in the following way:
% 
% signal = ([0 0, 0 0] - [0 0, 0 90]) - ([0 90, 0 0] + [0 90, 0 90]);
%                  
% where, in [X1 X2, Y1 Y2], X1 and X2 are the phases of the first and
% second refocusing pulses along X gradient, respectively, and Y1 and Y2
% are the phase of the third and fourth refocusing pulses along Y gradient.
% 
% Finally, this code simulates the spectrum at a given point in space (x,y),
% given the values of the slice selection gradients (Gx, and Gy).  The pulse
% waveform is assumed to be the same for both refocusing pulses.  In order
% to fully simulate the sLASER experiment, you have to run this
% simulation many times at various points in space (x,y), and then summate
% and scale together the resulting spectra.  
% 
% Feb 2020 - Jamie Near:  This code now accepts gradient modulated pulses.  
% If the input pulse is gradient modulated (waveform has 4 columns), then 
% the input parameters Gx and Gy are scaling factors for the GM function
% in order to achieve the desired slice thickness.
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% te        = echo time of sLASER experiment (ms)
% RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = If RF is not a gradient modulated pulse, Gx is the gradient strength
%             for first selective refocusing pulse [G/cm].  If RF is a gradient
%             modulated pulse, then Gx is a scaling factor for the GM
%             function to achieve the desired slice thickness in the
%             x-direction.  
% Gy        = If RF is not a gradient modulated pulse, Gy is the gradient strength
%             for first selective refocusing pulse [G/cm].  If RF is a gradient
%             modulated pulse, then Gy is a scaling factor for the GM
%             function to achieve the desired slice thickness in the
%             y-direction.  
% phl    = initial phase of the first refocusing pulse in [degrees];
% ph2   = initial phase of the second refocusing pulse in [degrees];
% ph3   = initial phase of the third refocusing pulse in [degrees];
% ph4   = initial phase of the fourth refocusing pulse in [degrees];
% flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
% centreFreq= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.

function out = sim_semiLASER_shaped_phCyc(n,sw,Bfield,linewidth,sys,te,RF,tp,dx,dy,Gx,Gy,ph1,ph2,ph3,ph4,flipAngle,centreFreq)

if nargin<18
    centreFreq=2.3;
    if nargin<17
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, scale the GM functions
% according to Gx and Gy, and then set Gx and Gy both equal to zero:
if RF.isGM
    RFX=rf_scaleGrad(RF,Gx);
    RFY=rf_scaleGrad(RF,Gy);
    Gx=0;
    Gy=0;
else
    RFX=RF;
    RFY=RF;
end

 if (te/4)<(tp/1000)
     error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
 end

%initialize evolution times
 tau1=(te/4-tp)/2;
 tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************ 
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
d=sim_shapedRF(d,H,RFX,tp,flipAngle,ph1,dx,Gx);         %1st shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RFX,tp,flipAngle,ph2,dx,Gx);         %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RFY,tp,flipAngle,ph3,dy,Gy);         %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RFY,tp,flipAngle,ph4,dy,Gy);         %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1

[out,dout]=sim_readout(d,H,n,sw,linewidth,90);          %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='semi-LASER';
out.te=te;
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


end
