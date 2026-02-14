% rf2d_getRFPulseEnergy.m
% Jamie Near, 2025.
%
% USAGE:
%   E = rf2d_getRFPulseEnergy(rf2d, peakPower)
%
% DESCRIPTION:
%   Compute the total energy [kJ] delivered by a 2D RF pulse by integrating
%   the instantaneous transmit power over the pulse duration. The power
%   waveform is taken proportional to |B1(t)|^2 and normalized to the peak
%   of the provided B1 waveform. If the instantaneous peak transmit power
%   at max |B1| is unknown, estimate it with rf2d_getRFPeakPower.
%   Useful for comparing total energy (and SAR proxies for equal durations)
%   across different RF pulse waveforms.
%
% INPUTS:
%   rf2d      - 2D RF pulse structure with field:
%                 • waveform : B1(t) samples (real or complex).
%                 • tp       : total pulse duration [s].    
%   peakPower - Instantaneous peak transmit power [kW] at max |B1|.
%
% OUTPUTS:
%   E         - Total pulse energy [kJ].
%
% NOTES:
%   • The waveform is internally normalized by max |B1|.
%   • See also: rf2d_getRFPeakPower


function [E]=rf2d_getRFPulseEnergy(rf2d,peakPower)

%Calculate the power waveform as the square of the normalized B1 waveform,
%multiplied by the peak power of the pulse.
powerWaveform = peakPower*(rf2d.waveform.^2)/max(rf2d.waveform.^2);

%calculate the step size of the RF pulse
deltaT=rf2d.tp/length(rf2d.waveform);

%Calculate the total pulse energy as the sum of the power waveform,
%multiplied by the step size:
E=sum(powerWaveform)*deltaT;

