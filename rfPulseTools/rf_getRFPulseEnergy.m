% rf_getRFPulseEnergy.m
% Jamie Near, McGill University 2021.
%
% USAGE:
% [E]=rf_getRFPulseEnergy(RF,tp,peakPower);
% 
% DESCRIPTION:
% Calculate the total pulse energy for an RF pulse, given its waveform and  
% duration, and given the instantaneous power delivered during the 
% peak amplitude of the pulse waveform.  If you do not know the instantaneous 
% peak power of the RF pulse, you can estimate it using rf_getRFPeakPower,
% which calculates the peak power based on the required B1 of the pulse,
% the peak B1 of the MRI system, and the system Tx amplifier power rating.
% This function is useful for comparing the specific absorption rate of 
% different RF pulse waveforms.
%
% 
% INPUTS:
% RF        = RF pulse definition structure
% tp        = pulse duration in [ms]
% peakPower = Instantaneous peak power [kW] of the RF pulse.  If this is
%             not known, you can use rf_getRFPeakPower to estimate.  
%
% OUTPUTS:
% E         = Total pulse energy in [kJ]


function [E]=rf_getRFPulseEnergy(RF,tp,peakPower);

%Calculate the power waveform as the square of the normalized B1 waveform,
%multiplied by the peak power of the pulse.
powerWaveform = peakPower*(RF.waveform.^2)/max(RF.waveform.^2);

%calculate the step size of the RF pulse
deltaT=tp/1000/length(RF.waveform);

%Calculate the total pulse energy as the sum of the power waveform,
%multiplied by the step size:
E=sum(powerWaveform)*deltaT;

