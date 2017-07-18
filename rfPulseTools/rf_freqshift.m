%rf_freqshift.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF_shift=rf_freqshift(RF,Tp,F);
% 
% DESCRIPTION:
% Apply a frequency shift to an RF pulse.
% 
% INPUTS:
% RF         = RF pulse definition structure.
% Tp         = duration of the rf pulse in [ms].
% F          = amount that you would like to frequency shift the rf pulse in [Hz].
%
% OUTPUTS:
% RF_shift   = Output rf pulse following frequency shift.

function RF_shift=rf_freqshift(RF,Tp,F)

N=size(RF.waveform,1);
Tp=Tp/1000;
dt=Tp/N;
t=[0:dt:Tp-dt];

phaseRamp = t * F * 360;

RF_shift=RF;
RF_shift.waveform(:,1)=RF_shift.waveform(:,1)+phaseRamp';
RF_shift.f0=RF.f0+F;

